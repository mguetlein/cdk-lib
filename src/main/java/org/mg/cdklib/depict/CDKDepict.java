package org.mg.cdklib.depict;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.SwingConstants;

import org.apache.batik.swing.JSVGCanvas;
import org.apache.commons.io.FileUtils;
import org.mg.cdklib.AtomContainerUtil;
import org.mg.cdklib.CDKConverter;
import org.mg.imagelib.ImageLoader;
import org.mg.imagelib.ImageType;
import org.mg.imagelib.ImageUtil;
import org.mg.imagelib.MultiImageIcon;
import org.mg.imagelib.MultiImageIcon.Layout;
import org.mg.imagelib.MultiImageIcon.Orientation;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.CollectionUtil;
import org.mg.javalib.util.SwingUtil;
import org.openscience.cdk.depict.Depiction;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

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

	@SuppressWarnings("unchecked")
	private static BufferedImage depictConnected(IAtomContainer mol) throws CDKException
	{
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(mol);
		sdg.generateCoordinates();
		mol = sdg.getMolecule();

		@SuppressWarnings("rawtypes")
		List generators = new ArrayList<>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		generators.add(new BasicAtomGenerator());
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators,
				new AWTFontManager());

		renderer.setup(mol, new Rectangle(0, 0, 1, 1));
		Rectangle diagramRectangle = renderer.paint(mol, new AWTDrawVisitor(
				(Graphics2D) new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB).getGraphics()));
		int width = (int) (1.5 * (diagramRectangle.getWidth() + diagramRectangle.x));
		int height = (int) (1.5 * (diagramRectangle.getHeight() + diagramRectangle.y));
		width = Math.max(20, width);
		height = Math.max(20, height);

		Rectangle drawArea = new Rectangle(width, height);
		renderer.setup(mol, drawArea);
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2 = (Graphics2D) image.getGraphics();
		//		g2.setColor(Color.WHITE);
		//		g2.fillRect(0, 0, width, height);
		renderer.paint(mol, new AWTDrawVisitor(g2), drawArea, true);

		return image;
	}

	public static void depictMatchToPNG(String pngFile, IAtomContainer mol, Color[] cols,
			boolean crop, int size) throws Exception
	{
		BufferedImage image = depictMatch(mol, cols, crop, size);
		ImageIO.write((RenderedImage) image, "PNG", new File(pngFile));
	}

	public static void depictMatchToPNG(String pngFile, IAtomContainer mol, int atoms[],
			boolean highlightOutgoingBonds, Color col, boolean crop, int maxSize) throws Exception
	{
		BufferedImage image = depictMatch(mol, atoms, highlightOutgoingBonds, col, crop, maxSize);
		ImageIO.write((RenderedImage) image, "PNG", new File(pngFile));
	}

	public static String COLOR_PROP = "colorProp";

	public static BufferedImage depictMatch(IAtomContainer mol, int atoms[],
			boolean highlightOutgoingBonds, Color col, boolean crop, int size)
			throws CloneNotSupportedException, CDKException
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
		return depictMatch(mol, new Color[] { col }, crop, size);
	}

	public static BufferedImage depictMatch(IAtomContainer mol, Color palette[], boolean crop,
			int size) throws CDKException
	{
		BufferedImage image = null;
		IAtomContainerSet set = ConnectivityChecker.partitionIntoMolecules(mol);
		if (set.getAtomContainerCount() < 2)
			image = depictMatchConnected(mol, palette, crop, size);
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
					image = depictMatchConnected(set.getAtomContainer(i), palette, true, size);
					break;
				}
			}
		}
		else
		{
			List<ImageIcon> icons = new ArrayList<ImageIcon>();
			for (int i = 0; i < set.getAtomContainerCount(); i++)
				icons.add(new ImageIcon(depictMatchConnected(set.getAtomContainer(i), palette, crop,
						crop ? size : -1)));
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
			boolean crop, int size) throws CDKException
	{
		if (crop && size == -1)
			throw new IllegalArgumentException();

		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(mol);
		sdg.generateCoordinates();
		mol = sdg.getMolecule();
		HashMap<IChemObject, Integer> ids = new HashMap<IChemObject, Integer>();
		IAtomContainer s = new AtomContainer();
		for (int j = 0; j < mol.getAtomCount(); j++)
			if (mol.getAtom(j).getProperty(COLOR_PROP, Integer.class) != null)
			{
				ids.put(mol.getAtom(j), mol.getAtom(j).getProperty(COLOR_PROP, Integer.class) + 1);
				s.addAtom(mol.getAtom(j));
			}
			else
				ids.put(mol.getAtom(j), 0);
		for (int j = 0; j < mol.getBondCount(); j++)
			if (mol.getBond(j).getProperty(COLOR_PROP, Integer.class) != null)
				ids.put(mol.getBond(j), mol.getBond(j).getProperty(COLOR_PROP, Integer.class) + 1);

		mol.setProperty(HighlightGenerator.ID_MAP, ids);

		@SuppressWarnings("rawtypes")
		List generators = new ArrayList<>();
		generators.add(new BasicSceneGenerator());
		generators.add(new HighlightGenerator());
		generators.add(new BasicBondGenerator());
		generators.add(new BasicAtomGenerator());
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators,
				new AWTFontManager());

		renderer.getRenderer2DModel().set(HighlightGenerator.HighlightPalette.class,
				HighlightGenerator.createPalette(null, palette));
		renderer.getRenderer2DModel().set(HighlightGenerator.HighlightRadius.class, 12.5);

		// determine preferred size of molecule
		renderer.setup(mol, new Rectangle(0, 0, 1, 1));
		Rectangle diagramRectangle = renderer.paint(mol, new AWTDrawVisitor(
				(Graphics2D) new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB).getGraphics()));
		int width = (int) (1.5 * (diagramRectangle.getWidth() + diagramRectangle.x));
		int height = (int) (1.5 * (diagramRectangle.getHeight() + diagramRectangle.y));
		width = Math.max(20, width);
		height = Math.max(20, height);

		// draw according to preferred size (with 10 pixels extra for the highlights)
		width += 10;
		height += 10;
		Rectangle drawArea = new Rectangle(5, 5, width - 10, height - 10);
		renderer.setup(mol, drawArea);
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2 = (Graphics2D) image.getGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, width, height);
		renderer.paint(mol, new AWTDrawVisitor(g2), drawArea, true);

		final Rectangle2D r = renderer.calculateDiagramBounds(s);
		int x, y, w, h;
		if (size >= width)
		{
			x = 0;
			w = width - 1;
		}
		else
		{
			x = (int) (r.getX() - (size - r.getWidth()) / 2.0);
			if (x + size > width)
				x -= (x + size - width);
			if (x < 0)
				x = 0;
			w = size;
		}
		if (size >= height)
		{
			y = 0;
			h = height - 1;
		}
		else
		{
			y = (int) (r.getY() - (size - r.getHeight()) / 2.0);
			if (y + size > height)
				y -= (y + size - height);
			if (y < 0)
				y = 0;
			h = size;
		}
		final Rectangle2D r2 = r.getFrame();
		r2.setFrame(x, y, w, h);
		//		System.out.println(r);
		//		System.out.println(r2);
		//		System.out.println(image.getWidth() + ", " + image.getHeight());

		//		JLabel l = new JLabel(new ImageIcon(image))
		//		{
		//			@Override
		//			public void paint(Graphics g)
		//			{
		//				super.paint(g);
		//				g.setColor(Color.GREEN);
		//				g.drawRect(0, 0, width, height);
		//				g.setColor(Color.RED);
		//				g.drawRect((int) r.getX(), (int) r.getY(), (int) r.getWidth(), (int) r.getHeight());
		//				g.setColor(Color.BLUE);
		//				g.drawRect((int) r2.getX(), (int) r2.getY(), (int) r2.getWidth(), (int) r2.getHeight());
		//
		//			}
		//		};
		//		SwingUtil.showInDialog(l);

		if (crop)
			image = image.getSubimage(x, y, w, h);
		if (size != -1 && (image.getWidth() > size || image.getHeight() > size))
			image = (BufferedImage) ImageLoader.getShrinkedImage(new ImageIcon(image), size, size)
					.getImage();
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
		return depictMatch(mol, a, false, col, crop, size);
	}

	public static byte[] depictMatchToPNGBytes(IAtomContainer mol, String smarts, Color col,
			boolean crop, int size) throws IOException, CDKException, CloneNotSupportedException
	{
		return ImageUtil.imageToPNGBytes(depictMatch(mol, smarts, col, crop, size));
	}

	private static Depiction depictMatch2Depiction(IAtomContainer mol, String smarts)
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

	public static BufferedImage depictMatch2Image(IAtomContainer mol, String smarts)
			throws IOException, CDKException
	{
		return depictMatch2Depiction(mol, smarts).toImg();
	}

	public static byte[] depictMatch2Bytes(IAtomContainer mol, String smarts, ImageType type)
			throws IOException, CDKException
	{
		if (type == ImageType.svg)
			return depictMatch2Depiction(mol, smarts).toSvgStr().getBytes();
		else
			return ImageUtil.imageToPNGBytes(depictMatch2Image(mol, smarts));
	}

	public static void depictMatch2File(IAtomContainer mol, String smarts, String file)
			throws IOException, CDKException, CloneNotSupportedException
	{
		FileUtils.writeByteArrayToFile(new File(file),
				depictMatch2Bytes(mol, smarts, ImageType.fromFilename(file)));
	}

	public static JSVGCanvas depictMatch2ToSVGCanvas(IAtomContainer mol, String smarts)
			throws IOException, CDKException, CloneNotSupportedException
	{
		File f = File.createTempFile("smiles", ".svg");
		depictMatch2File(mol, smarts, f.getPath());
		return ImageUtil.svgFileToCanvas(f.getPath());
	}

	public static void demo() throws InvalidSmilesException, Exception
	{

		//		SwingUtil.showInDialog(
		//				new JLabel(new ImageIcon(depictMatch2(CDKConverter.parseSmiles("CCC=O"), null))));
		SwingUtil.showInDialog(depictMatch2ToSVGCanvas(
				CDKConverter.parseSmiles(
						"OC1C2C(N(C)C)C(=O)C(=C(O)N)C(=O)C2(O)C(=O)C2=C(O)c3c(C(C12)(C)O)c(Cl)ccc3O"),
				null));
		//		{
		//			DefaultFormBuilder b = new DefaultFormBuilder(new FormLayout("p"));
		//			for (Integer size : new Integer[] { 100, -1 })
		//			{
		//				for (String smiles : new String[] { "[Na+].[Na+].O=S(C1=CC=C(C(C)CC)C=C1)([O-])=O",
		//						"CC(C)C(C1=CC=C(C=C1)Cl)C(=O)OC(C#N)C2=CC(=CC=C2)OC3=CC=CC=C3" })
		//				{
		//					BufferedImage img = draw(
		//							new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles), size);
		//					b.append(getLabel(img, smiles + " size:" + size));
		//				}
		//			}
		//			SwingUtil.showInFrame(new JScrollPane(b.getPanel()));
		//		}
		{
			//			DefaultFormBuilder b = new DefaultFormBuilder(new FormLayout("p,20dlu,p,20dlu,p"));
			//			int atoms[] = new int[] { 2, 3, 4 };
			//			for (Integer size : new Integer[] { 75, -1 })
			//			{
			//				for (boolean crop : new boolean[] { true, false })
			//				{
			//					if (crop && size == -1)
			//						continue;
			//
			//					for (boolean bonds : new boolean[] { true, false })
			//					{
			//						for (String smiles : new String[] { "O=C1C2=C(C=CC=C2)C(=O)C3=C1C=CC=C3",
			//								"C1(=C(C=CC=C1)N)OC.[H]Cl", "[Na].[H]Cl[Mg]" })
			//						{
			//
			//							BufferedImage img = depictMatch(
			//									new SmilesParser(SilentChemObjectBuilder.getInstance())
			//											.parseSmiles(smiles),
			//									atoms, bonds, Color.RED, crop, size);
			//							b.append(getLabel(img, smiles + " crop:" + crop + " size:" + size
			//									+ " bonds:" + bonds));
			//						}
			//					}
			//				}
			//			}
			//			SwingUtil.showInFrame(new JScrollPane(b.getPanel()));
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
	}

}
