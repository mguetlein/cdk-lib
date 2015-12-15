package org.mg.cdklib;

import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

public class CDKTesting
{
	public static void main(String[] args) throws CDKException
	{
		{
			//String smi2 = "Cl";
			Graphics.draw("Cl", "[H]Cl");
		}

		String smi = "CCC(=O)O";
		IAtomContainer mol = CDKConverter.parseSmiles(smi);
		System.out.println(new SmilesGenerator().create(mol));
		dehydrolyse(mol);
		System.out.println(new SmilesGenerator().create(mol));

		//		String smi = "CCC(=O)O";
		//		IAtomContainer mol = Converter.parseSmiles(smi);
		//		String carboxylSMARTS = "C(O)=O";
		//		int oIndexInSMARTS = 1;
		//		SMARTSQueryTool smartsMatch = new SMARTSQueryTool(carboxylSMARTS, SilentChemObjectBuilder.getInstance());
		//		if (smartsMatch.matches(mol))
		//		{
		//			for (List<Integer> match : smartsMatch.getMatchingAtoms())
		//			{
		//				IAtom oAtom = mol.getAtom(match.get(oIndexInSMARTS));
		//				oAtom.setFormalCharge(-1);
		//				oAtom.setImplicitHydrogenCount(0);
		//			}
		//		}
		//		System.out.println(new SmilesGenerator().create(mol));
		//
		//		Graphics.draw(new SmilesGenerator().create(mol));
	}

	public static void dehydrolyse(IAtomContainer mol)
	{
		HashMap<String, Integer> smartsToOxIdx = new HashMap<>();
		smartsToOxIdx.put("C(O)=O", 1);
		for (String smarts : smartsToOxIdx.keySet())
		{
			try
			{
				SMARTSQueryTool smartsMatch = new SMARTSQueryTool(smarts, SilentChemObjectBuilder.getInstance());
				if (smartsMatch.matches(mol))
				{
					for (List<Integer> match : smartsMatch.getMatchingAtoms())
					{
						IAtom oxygen = mol.getAtom(match.get(smartsToOxIdx.get(smarts)));
						oxygen.setFormalCharge(-1);
						oxygen.setImplicitHydrogenCount(0);
					}
				}
			}
			catch (CDKException e)
			{
				throw new RuntimeException(e);
			}
		}
	}
}
