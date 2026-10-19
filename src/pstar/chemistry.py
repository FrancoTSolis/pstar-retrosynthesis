"""Small RDKit utilities shared by policies and validators."""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def canonicalize(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


@lru_cache(maxsize=200_000)
def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol is not None else 50


@lru_cache(maxsize=200_000)
def basic_descriptors(smiles: str) -> tuple[float, ...]:
    """Six descriptors that do not require external learned models."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (5.0, 3.0, 2.0, 2.0, 2.0, 2.0)
    return (
        mol.GetNumHeavyAtoms() / 10.0,
        float(rdMolDescriptors.CalcNumRings(mol)),
        float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        float(len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))),
        rdMolDescriptors.CalcNumHeteroatoms(mol) / 5.0,
        rdMolDescriptors.CalcNumRotatableBonds(mol) / 5.0,
    )
