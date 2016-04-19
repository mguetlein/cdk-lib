package org.mg.cdklib;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.List;

import org.mg.javalib.io.External3DComputer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class CDKTesting
{
	public static void normalizeTest() throws CDKException
	{
		normalize("C(=O)O");
		normalize("C(=O)[O-]");
	}

	public static void normalize(String smi) throws CDKException
	{
		IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
				.parseSmiles(smi);

		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

		System.out.println(new SmilesGenerator().create(mol));

	}

	public static void main(String[] args) throws CDKException, IOException
	{
		normalizeTest();

		//		parseMol();

		//		IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
		//				.parseSmiles("COC1=C2C=CC=CC2=CC(=C1O)O");
		//		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
		//				Cycles.or(Cycles.all(), Cycles.edgeShort()));
		//		aromaticity.apply(mol);
		//		System.out.println(countAromAtoms(mol));

		//		{
		//			//String smi2 = "Cl";
		//			Graphics.draw("Cl", "[H]Cl");
		//		}
		//
		//		String smi = "CCC(=O)O";
		//		IAtomContainer mol = CDKConverter.parseSmiles(smi);
		//		System.out.println(new SmilesGenerator().create(mol));
		//		dehydrolyse(mol);
		//		System.out.println(new SmilesGenerator().create(mol));

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

	public static void parseMol() throws CDKException, IOException
	{
		String mol = External3DComputer.get3D("CCC");
		ISimpleChemObjectReader reader = new MDLV2000Reader(
				new ByteArrayInputStream(mol.getBytes(StandardCharsets.UTF_8)));
		IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
		IAtomContainer moleclue = (IAtomContainer) ChemFileManipulator.getAllAtomContainers(content)
				.get(0);
		reader.close();
		moleclue = (IAtomContainer) AtomContainerManipulator.removeHydrogens(moleclue);
		System.out.println(moleclue.getAtomCount());
		System.out.println(moleclue.getAtom(0).getPoint3d());
	}

	public static void isotope() throws CDKException
	{
		IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
				.parseSmiles("[14C]");
		System.out.println(CDKConverter.toSmiles(mol));

	}

	public static int countAromAtoms(IAtomContainer mol)
	{
		int c = 0;
		for (IAtom atom : mol.atoms())
			if (atom.getFlag(CDKConstants.ISAROMATIC))
				c++;
		return c;
	}

	public static void dehydrolyse(IAtomContainer mol)
	{
		HashMap<String, Integer> smartsToOxIdx = new HashMap<>();
		smartsToOxIdx.put("C(O)=O", 1);
		for (String smarts : smartsToOxIdx.keySet())
		{
			try
			{
				SMARTSQueryTool smartsMatch = new SMARTSQueryTool(smarts,
						SilentChemObjectBuilder.getInstance());
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
