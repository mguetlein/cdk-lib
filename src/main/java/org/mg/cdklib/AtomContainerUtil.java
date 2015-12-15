package org.mg.cdklib;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;

public class AtomContainerUtil
{
	public static List<IChemObject> getAtomsAndBonds(IAtomContainer mol)
	{
		List<IChemObject> l = new ArrayList<>();
		for (int i = 0; i < mol.getAtomCount(); i++)
			l.add(mol.getAtom(i));
		for (int i = 0; i < mol.getBondCount(); i++)
			l.add(mol.getBond(i));
		return l;
	}
}
