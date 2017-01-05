package org.mg.cdklib.depict;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.mg.cdklib.AtomContainerUtil;
import org.mg.imagelib.ImageLoader;
import org.mg.imagelib.ImageUtil;
import org.mg.imagelib.MultiImageIcon;
import org.mg.imagelib.MultiImageIcon.Layout;
import org.mg.imagelib.MultiImageIcon.Orientation;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.CollectionUtil;
import org.mg.javalib.util.SwingUtil;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.SymbolVisibility;
import org.openscience.cdk.renderer.color.CDK2DAtomColors;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.standard.SelectionVisibility;
import org.openscience.cdk.renderer.generators.standard.StandardGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.FormLayout;

/**
 * used by coffer and cfpminer
 * * custom splitting of unconnected molecules
 * * cropping and higlighting of selected atoms
 * * custom highlighter OuterGlowWhiteEdge 
 * 
 * @author martin
 *
 */
public class CDKDepict
{
	public static void depictToPNG(String pngFile, IAtomContainer mol, int maxSize)
			throws CDKException, IOException
	{
		BufferedImage image = depict(mol, maxSize);
		ImageIO.write((RenderedImage) image, "PNG", new File(pngFile));
	}

	public static byte[] depictToPNGBytes(IAtomContainer mol, int maxSize)
			throws CDKException, IOException
	{
		return ImageUtil.imageToPNGBytes(depict(mol, maxSize));
	}

	public static BufferedImage depict(IAtomContainer mol, int maxSize) throws CDKException
	{
		IAtomContainerSet set = ConnectivityChecker.partitionIntoMolecules(mol);
		BufferedImage image;
		if (set.getAtomContainerCount() < 2)
			image = depictConnected(mol);
		else
		{
			List<ImageIcon> icons = new ArrayList<ImageIcon>();
			for (int i = 0; i < set.getAtomContainerCount(); i++)
				icons.add(new ImageIcon(depictConnected(set.getAtomContainer(i))));
			image = (BufferedImage) new MultiImageIcon(icons, Layout.horizontal, Orientation.center,
					2).getImage();
		}
		if (maxSize != -1)
			image = (BufferedImage) ImageLoader
					.getShrinkedImage(new ImageIcon(image), maxSize, maxSize).getImage();
		return image;
	}

	private static BufferedImage depictConnected(IAtomContainer mol) throws CDKException
	{
		return new DepictionGenerator().withTerminalCarbons().withAtomColors().depict(mol).toImg();
	}

	public static void depictMatchToPNG(String pngFile, IAtomContainer mol, Color[] cols,
			boolean crop, int size, boolean drawCarbonSymbolsInSelection) throws Exception
	{
		BufferedImage image = depictMatch(mol, cols, crop, size, drawCarbonSymbolsInSelection);
		ImageIO.write((RenderedImage) image, "PNG", new File(pngFile));
	}

	public static void depictMatchToPNG(String pngFile, IAtomContainer mol, int atoms[],
			boolean highlightOutgoingBonds, Color col, boolean crop, int maxSize,
			boolean drawCarbonSymbolsInSelection) throws Exception
	{
		BufferedImage image = depictMatch(mol, atoms, highlightOutgoingBonds, col, crop, maxSize,
				drawCarbonSymbolsInSelection);
		ImageIO.write((RenderedImage) image, "PNG", new File(pngFile));
	}

	public static String COLOR_PROP = "colorProp";

	public static BufferedImage depictMatch(IAtomContainer mol, int atoms[],
			boolean highlightOutgoingBonds, Color col, boolean crop, int size,
			boolean drawCarbonSymbolsInSelection) throws CloneNotSupportedException, CDKException
	{
		if (atoms == null || atoms.length == 0)
			if (crop)
				throw new IllegalArgumentException();
		mol = mol.clone();
		for (IChemObject c : AtomContainerUtil.getAtomsAndBonds(mol))
			c.setProperty(COLOR_PROP, null);
		for (int j = 0; j < mol.getAtomCount(); j++)
		{
			if (atoms == null || ArrayUtil.indexOf(atoms, j) == -1)
				continue;
			mol.getAtom(j).setProperty(COLOR_PROP, 0);
			for (int i = 0; i < mol.getAtomCount(); i++)
			{
				if (i == j)
					continue;
				if (!highlightOutgoingBonds && ArrayUtil.indexOf(atoms, i) == -1)
					continue;
				IBond b = mol.getBond(mol.getAtom(j), mol.getAtom(i));
				if (b != null)
					b.setProperty(COLOR_PROP, 0);
			}
		}
		return depictMatch(mol, new Color[] { col }, crop, size, drawCarbonSymbolsInSelection);
	}

	public static BufferedImage depictMatch(IAtomContainer mol, Color palette[], boolean crop,
			int size, boolean drawCarbonSymbolsInSelection) throws CDKException
	{
		BufferedImage image = null;
		IAtomContainerSet set = ConnectivityChecker.partitionIntoMolecules(mol);
		if (set.getAtomContainerCount() < 2)
			image = depictMatchConnected(mol, palette, crop, size, drawCarbonSymbolsInSelection);
		else if (crop)
		{
			for (int i = 0; i < set.getAtomContainerCount(); i++)
			{
				boolean match = false;
				for (IChemObject c : AtomContainerUtil.getAtomsAndBonds(set.getAtomContainer(i)))
					if (c.getProperty(COLOR_PROP, Integer.class) != null)
						match = true;
				if (match)
				{
					image = depictMatchConnected(set.getAtomContainer(i), palette, true, size,
							drawCarbonSymbolsInSelection);
					break;
				}
			}
		}
		else
		{
			List<ImageIcon> icons = new ArrayList<ImageIcon>();
			for (int i = 0; i < set.getAtomContainerCount(); i++)
				icons.add(new ImageIcon(depictMatchConnected(set.getAtomContainer(i), palette, crop,
						crop ? size : -1, drawCarbonSymbolsInSelection)));
			image = (BufferedImage) new MultiImageIcon(icons, Layout.horizontal, Orientation.center,
					2).getImage();
		}
		if (!crop && size != -1)
			image = (BufferedImage) ImageLoader.getShrinkedImage(new ImageIcon(image), size, size)
					.getImage();
		return image;
	}

	/**
	 * atoms / bonds must have property {@link COLOR_PROP} assigned
	 * null -> no highlight
	 * value (must be integer) -> index in pallete color array
	 * 
	 * 
	 * @param mol
	 * @param palette
	 * @param crop
	 * @param size
	 * @return
	 * @throws CDKException 
	 */
	@SuppressWarnings("unchecked")
	private static BufferedImage depictMatchConnected(IAtomContainer mol, Color palette[],
			boolean crop, int size, boolean drawCarbonSymbolsInSelection) throws CDKException
	{
		if (crop && size == -1)
			throw new IllegalArgumentException();

		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(mol);
		sdg.generateCoordinates();
		mol = sdg.getMolecule();

		IAtomContainer matchingAtoms = new AtomContainer();
		// set color of atoms
		for (int j = 0; j < mol.getAtomCount(); j++)
		{
			IAtom a = mol.getAtom(j);
			if (a.getProperty(COLOR_PROP, Integer.class) != null)
			{
				// set property for renderer
				a.setProperty(StandardGenerator.HIGHLIGHT_COLOR,
						palette[a.getProperty(COLOR_PROP, Integer.class)]);
				matchingAtoms.addAtom(a);
			}
			else
			{
				// if no color is set, set color to white only if a connected bond is highlighted
				// this is done because the carbons of selected atoms are shown explicitly (if enabled)
				for (IBond b : mol.getConnectedBondsList(a))
				{
					if (b.getProperty(COLOR_PROP, Integer.class) != null)
					{
						a.setProperty(StandardGenerator.HIGHLIGHT_COLOR, Color.WHITE);
						break;
					}
				}
			}
		}
		for (int j = 0; j < mol.getBondCount(); j++)
		{
			IBond b = mol.getBond(j);
			if (b.getProperty(COLOR_PROP, Integer.class) != null)
				b.setProperty(StandardGenerator.HIGHLIGHT_COLOR,
						palette[b.getProperty(COLOR_PROP, Integer.class)]);
		}

		@SuppressWarnings("rawtypes")
		List generators = new ArrayList<>();
		generators.add(new BasicSceneGenerator());
		generators.add(new StandardGenerator(new Font(Font.SANS_SERIF, Font.PLAIN, 18)));
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators,
				new AWTFontManager());

		if (drawCarbonSymbolsInSelection)
			renderer.getRenderer2DModel().set(StandardGenerator.Visibility.class,
					SelectionVisibility.all(SymbolVisibility.iupacRecommendations()));
		else
			renderer.getRenderer2DModel().set(StandardGenerator.Visibility.class,
					SymbolVisibility.iupacRecommendations());

		renderer.getRenderer2DModel().set(StandardGenerator.AtomColor.class, new CDK2DAtomColors());
		renderer.getRenderer2DModel().set(StandardGenerator.Highlighting.class,
				StandardGenerator.HighlightStyle.OuterGlowWhiteEdge);
		renderer.getRenderer2DModel().set(StandardGenerator.OuterGlowWidth.class, 4.0);

		// this is not nice, but its working, better would be to get the size from the renderer
		BufferedImage testImg = new DepictionGenerator().withCarbonSymbols().depict(mol).toImg();
		int width = testImg.getWidth();
		int height = testImg.getHeight();
		// add 20 pixels for the highlights
		int extra = 20;
		width += extra;
		height += extra;

		// draw image with fixed size
		final Rectangle drawArea = new Rectangle(0, 0, width, height);
		renderer.setup(mol, drawArea);
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2 = (Graphics2D) image.getGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, width, height);
		renderer.paint(mol, new AWTDrawVisitor(g2), drawArea, true);

		if (crop)
		{
			// calculate bounds of selection
			Rectangle bounds = renderer.calculateDiagramBounds(matchingAtoms);

			// add extra for extra rendered hydrogens that are unfortunately not in bounds
			extra = 30;
			final Rectangle r = new Rectangle(bounds.x - extra / 2, bounds.y - extra / 2,
					bounds.width + extra, bounds.height + extra);

			// size is always given when cropping
			// extend the selection to be quadratic and of minimum-size: size
			int minSize = Math.max(size, Math.max(r.width, r.height));
			int addToWidth = r.width < minSize ? minSize - r.width : 0;
			int addToHeight = r.height < minSize ? minSize - r.height : 0;
			final Rectangle r2 = new Rectangle(r.x - addToWidth / 2, r.y - addToHeight / 2,
					r.width + addToWidth, r.height + addToHeight);

			// selection should not be larger than the image
			final Rectangle r3 = new Rectangle(Math.max(0, r2.x), Math.max(0, r2.y),
					Math.min(width - 1, r2.width), Math.min(height - 1, r2.height));
			// selection should not extend the image
			if (r3.x + r3.width > width)
				r3.x = width - r3.width;
			if (r3.y + r3.height > height)
				r3.y = height - r3.height;

			//			JLabel l = new JLabel(new ImageIcon(image))
			//			{
			//				@Override
			//				public void paint(Graphics g)
			//				{
			//					super.paint(g);
			//					g.setColor(Color.GREEN);
			//					g.drawRect(0, 0, (int) drawArea.getWidth(), (int) drawArea.getHeight());
			//					g.setColor(Color.RED);
			//					g.drawRect(r.x, r.y, r.width, r.height);
			//					g.setColor(Color.BLUE);
			//					g.drawRect(r2.x, r2.y, r2.width, r2.height);
			//					g.setColor(Color.MAGENTA);
			//					g.drawRect(r3.x, r3.y, r3.width, r3.height);
			//
			//				}
			//			};
			//			SwingUtil.showInDialog(l);
			image = image.getSubimage(r3.x, r3.y, r3.width, r3.height);
		}
		if (size != -1 && (image.getWidth() > size || image.getHeight() > size))
		{
			image = (BufferedImage) ImageLoader.getShrinkedImage(new ImageIcon(image), size, size)
					.getImage();
		}
		if (crop && (image.getWidth() < size || image.getHeight() < size))
		{
			BufferedImage img = new BufferedImage(size, size, BufferedImage.TYPE_INT_ARGB);
			Graphics g = img.getGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, size, size);
			g.drawImage(image, (size - image.getWidth()) / 2, (size - image.getHeight()) / 2, null);
			image = img;
		}
		return image;
	}

	private static JLabel getLabel(BufferedImage img, String text)
	{
		BufferedImage imgW = new BufferedImage(img.getWidth(), img.getHeight(),
				BufferedImage.TYPE_INT_RGB);
		imgW.getGraphics().setColor(Color.WHITE);
		imgW.getGraphics().fillRect(0, 0, imgW.getWidth(), imgW.getHeight());
		imgW.getGraphics().drawImage(img, 0, 0, null);
		JLabel l = new JLabel(text, new ImageIcon(imgW), SwingConstants.CENTER);
		l.setHorizontalTextPosition(JLabel.CENTER);
		l.setVerticalTextPosition(JLabel.BOTTOM);
		return l;
	}

	public static BufferedImage depictMatch(IAtomContainer mol, String smarts, Color col,
			boolean crop, int size) throws CDKException, IOException, CloneNotSupportedException
	{
		SmartsPattern sp = SmartsPattern.create(smarts, SilentChemObjectBuilder.getInstance());
		Set<Integer> atoms = new HashSet<>();
		Mappings map = sp.matchAll(mol);
		for (int[] is : map)
			for (int i : is)
				atoms.add(i);
		int a[] = new int[0];
		if (atoms.size() > 0)
			a = ArrayUtil.toPrimitiveIntArray(CollectionUtil.toArray(atoms));
		return depictMatch(mol, a, false, col, crop, size, true);
	}

	public static byte[] depictMatchToPNGBytes(IAtomContainer mol, String smarts, Color col,
			boolean crop, int size) throws IOException, CDKException, CloneNotSupportedException
	{
		return ImageUtil.imageToPNGBytes(depictMatch(mol, smarts, col, crop, size));
	}

	public static void demo() throws InvalidSmilesException, Exception
	{
		{
			DefaultFormBuilder b = new DefaultFormBuilder(new FormLayout("p"));
			for (Integer size : new Integer[] { 100, -1 })
			{
				for (String smiles : new String[] { "[Na+].[Na+].O=S(C1=CC=C(C(C)CC)C=C1)([O-])=O",
						"CC(C)C(C1=CC=C(C=C1)Cl)C(=O)OC(C#N)C2=CC(=CC=C2)OC3=CC=CC=C3",
						"O=C1C2=C(C=CC=C2)C(=O)C3=C1C=CC=C3", "[H]Cl[Mg]" })
				{
					BufferedImage img = depict(
							new SmilesParser(SilentChemObjectBuilder.getInstance())
									.parseSmiles(smiles),
							size);
					b.append(getLabel(img, smiles + " size:" + size));
				}
			}
			SwingUtil.showInFrame(new JScrollPane(b.getPanel()), "no match", false);
		}
		{
			DefaultFormBuilder b = new DefaultFormBuilder(new FormLayout("p,20dlu,p,20dlu,p"));
			int atoms[] = new int[] { 2, 3, 4, 5, 6 };
			for (boolean crop : new boolean[] { // 
					true, //
					false //
			})
			{
				for (Integer size : new Integer[] { // 
						150, //
						-1 //
				})
				{
					if (crop && size == -1)
						continue;

					for (boolean bonds : new boolean[] { // 
							true, //
							false //
					})
					{
						for (String smiles : new String[] { // 
								"O=C1C2=C(C=CC=C2)C(=O)C3=C1C=CC=C3", //
								"C1(=C(C=CC=C1)N)OC.[H]Cl", // 
								"[Na].[H]Cl[Mg]", //
								//"[H]Cl[Mg]" //
						})
						{

							BufferedImage img = depictMatch(
									new SmilesParser(SilentChemObjectBuilder.getInstance())
											.parseSmiles(smiles),
									atoms, bonds, Color.RED, crop, size, true);
							b.append(getLabel(img, "<html>" + smiles + "<br>crop:" + crop + " size:"
									+ size + " bonds:" + bonds + "</html>"));
						}
					}
				}
			}
			SwingUtil.showInFrame(new JScrollPane(b.getPanel()));
		}
		System.exit(1);

	}

	public static void main(String[] args) throws InvalidSmilesException, Exception
	{
		//		IAtomContainer mol = CDKConverter.parseSmiles("O");
		//		JPanel p = new JPanel();
		//		p.add(new JLabel(new ImageIcon(depict(mol, 150))));
		//		p.add(new JLabel(
		//				new ImageIcon(depictMatch(mol, new int[] { 0 }, false, Color.RED, true, 150))));
		//		p.add(new JLabel(
		//				new ImageIcon(depictMatch(mol, new int[] { 0 }, false, Color.RED, false, 150))));
		//		SwingUtil.showInFrame(p);
		//		SwingUtil.waitWhileWindowsVisible();
		//		System.exit(1);

		demo();
		System.exit(1);
	}

}
