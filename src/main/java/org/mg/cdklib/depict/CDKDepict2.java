package org.mg.cdklib.depict;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import org.apache.batik.swing.JSVGCanvas;
import org.apache.commons.io.FileUtils;
import org.mg.cdklib.CDKConverter;
import org.mg.imagelib.ImageType;
import org.mg.imagelib.ImageUtil;
import org.mg.javalib.util.SwingUtil;
import org.openscience.cdk.depict.Depiction;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * uses DepictionGenerator (instead of custom stuff)
 * 
 * @author martin
 *
 */
public class CDKDepict2
{
	static Depiction depictMatchToDepiction(IAtomContainer mol, String smarts)
			throws IOException, CDKException
	{
		AtomContainerManipulator.suppressHydrogens(mol);
		DepictionGenerator dptgen = new DepictionGenerator();
		if (smarts != null)
		{
			Pattern ptrn = SmartsPattern.create(smarts, SilentChemObjectBuilder.getInstance());
			Iterable<IChemObject> hits = ptrn.matchAll(mol).uniqueAtoms().toChemObjects();
			dptgen = dptgen.withHighlight(hits, Color.RED);
		}
		else
		{
			dptgen = dptgen.withAtomColors();
		}
		return dptgen.depict(mol);
	}

	public static BufferedImage depictMatchToImage(IAtomContainer mol, String smarts)
			throws IOException, CDKException
	{
		return depictMatchToDepiction(mol, smarts).toImg();
	}

	public static byte[] depictMatchToBytes(IAtomContainer mol, String smarts, ImageType type)
			throws IOException, CDKException
	{
		if (type == ImageType.svg)
			return depictMatchToDepiction(mol, smarts).toSvgStr().getBytes();
		else
			return ImageUtil.imageToPNGBytes(depictMatchToImage(mol, smarts));
	}

	public static void depictMatchToFile(IAtomContainer mol, String smarts, String file)
			throws IOException, CDKException, CloneNotSupportedException
	{
		FileUtils.writeByteArrayToFile(new File(file),
				depictMatchToBytes(mol, smarts, ImageType.fromFilename(file)));
	}

	public static JSVGCanvas depictMatchToSVGCanvas(IAtomContainer mol, String smarts)
			throws IOException, CDKException, CloneNotSupportedException
	{
		File f = File.createTempFile("smiles", ".svg");
		depictMatchToFile(mol, smarts, f.getPath());
		return ImageUtil.svgFileToCanvas(f.getPath());
	}

	public static void main(String[] args)
			throws InvalidSmilesException, IOException, CDKException, CloneNotSupportedException
	{
		//		SwingUtil.showInDialog(
		//				new JLabel(new ImageIcon(depictMatch2(CDKConverter.parseSmiles("CCC=O"), null))));
		SwingUtil.showInDialog(CDKDepict2.depictMatchToSVGCanvas(
				CDKConverter.parseSmiles(
						"OC1C2C(N(C)C)C(=O)C(=C(O)N)C(=O)C2(O)C(=O)C2=C(O)c3c(C(C12)(C)O)c(Cl)ccc3O"),
				null));
	}

}
