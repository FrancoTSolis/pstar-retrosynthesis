"""Adapter for the public DESP one-step template model."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from .chemistry import canonicalize


def _dependencies():
    try:
        import torch
        import torch.nn as nn
        from rdchiral.initialization import rdchiralReactants, rdchiralReaction
        from rdchiral.main import rdchiralRun
    except ImportError as error:
        raise ImportError('install the DESP extra with pip install -e ".[desp]"') from error
    return torch, nn, rdchiralReactants, rdchiralReaction, rdchiralRun


class DESPPredictor:
    def __init__(self, model_path: str | Path, templates_path: str | Path, device: str = "cpu"):
        torch, nn, _, _, _ = _dependencies()
        self.torch = torch
        with open(templates_path) as handle:
            self.templates = {int(key): value for key, value in json.load(handle).items()}
        checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
        args = checkpoint["args"]
        hidden = [int(value) for value in args.hidden_sizes.split(",")]
        activations = {
            "relu": nn.ReLU,
            "elu": nn.ELU,
            "gelu": nn.GELU,
            "leakyrelu": nn.LeakyReLU,
            "sigmoid": nn.Sigmoid,
            "tanh": nn.Tanh,
        }
        layers = []
        dimensions = [args.fp_size, *hidden]
        for source, target in zip(dimensions[:-1], dimensions[1:]):
            layers.extend([nn.Linear(source, target), activations[args.hidden_activation]()])
        self.model = nn.Sequential(*layers, nn.Linear(hidden[-1], args.n_templates))
        state = {key.replace("module.", ""): value for key, value in checkpoint["state_dict"].items()}
        translated = {}
        for key, value in state.items():
            if key.startswith("layers.") and ".linear." in key:
                parts = key.split(".")
                translated[f"{2 * int(parts[1])}.{parts[-1]}"] = value
            elif key.startswith("output_layer."):
                translated[f"{2 * len(hidden)}.{key.split('.')[-1]}"] = value
        self.model.load_state_dict(translated)
        self.model.to(device).eval()
        self.device = device
        self.fp_size = int(args.fp_size)
        self._compiled = {}

    def _fingerprint(self, smiles: str):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"invalid SMILES: {smiles}")
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(
            molecule,
            radius=2,
            nBits=self.fp_size,
            useChirality=True,
        )
        return np.asarray(fingerprint, dtype=np.float32)

    def predict(self, product: str, top_n: int = 50) -> list[dict]:
        torch, _, rdchiral_reactants, rdchiral_reaction, rdchiral_run = _dependencies()
        fingerprint = torch.from_numpy(self._fingerprint(product)).unsqueeze(0).to(self.device)
        with torch.no_grad():
            probabilities = torch.softmax(self.model(fingerprint), dim=1)
            scores, indices = torch.topk(probabilities, top_n)
        product_object = rdchiral_reactants(product)
        merged: dict[tuple[str, ...], float] = defaultdict(float)
        templates: dict[tuple[str, ...], str] = {}
        for index, score in zip(indices[0].tolist(), scores[0].tolist()):
            template = self.templates[index]
            reaction = self._compiled.get(template)
            if reaction is None:
                reaction_text = template
                if len(template.split(">>")[0].split(".")) > 1:
                    reaction_text = "(" + template.replace(">>", ")>>")
                reaction = rdchiral_reaction(reaction_text)
                self._compiled[template] = reaction
            try:
                outputs = rdchiral_run(reaction, product_object)
            except Exception:
                continue
            for output in outputs:
                reactants = tuple(
                    sorted(
                        canonical
                        for row in output.split(".")
                        if (canonical := canonicalize(row))
                    )
                )
                if not reactants:
                    continue
                merged[reactants] += float(score)
                templates.setdefault(reactants, template)
        total = sum(merged.values())
        return [
            {
                "reactants": list(reactants),
                "score": score / total if total else 0.0,
                "template": templates[reactants],
            }
            for reactants, score in sorted(merged.items(), key=lambda row: -row[1])
        ]


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate DESP one-step candidates")
    parser.add_argument("--model", required=True)
    parser.add_argument("--templates", required=True)
    parser.add_argument("--molecules", required=True, help="CSV with a smiles column")
    parser.add_argument("--output", required=True)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    with open(args.molecules, newline="") as handle:
        molecules = [row["smiles"] for row in csv.DictReader(handle)]
    predictor = DESPPredictor(args.model, args.templates, args.device)
    output = {}
    for smiles in molecules:
        canonical = canonicalize(smiles) or smiles
        output[canonical] = predictor.predict(canonical, args.top_n)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(output))
