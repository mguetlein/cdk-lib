package org.mg.cdklib;

import java.awt.Font;
import java.awt.GridLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.apache.commons.io.IOUtils;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.ImageLoader;
import org.mg.javalib.util.StringUtil;
import org.mg.javalib.util.SwingUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

public class Graphics
{
	public static final String CACHE_DIR = System.getProperty("user.home") + "/.cdk/cache/";

	public static String drawToFileOB(String smiles)
	{
		String type = "png";
		String file = CACHE_DIR + StringUtil.getMD5(smiles) + "." + type;
		if (!new File(file).exists())
		{
			try
			{
				String cmd[] = new String[] { "obabel", "-:" + smiles, "-o" + type };
				System.err.println("create img> " + ArrayUtil.toString(cmd, " "));
				System.err.flush();
				Process p = Runtime.getRuntime().exec(cmd);
				FileOutputStream fo = new FileOutputStream(new File(file));
				IOUtils.copy(p.getInputStream(), fo);
				fo.close();
				p.getInputStream().close();
				BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String s;
				while ((s = stdError.readLine()) != null)
					System.err.println("create img> " + s);
				stdError.close();
			}
			catch (IOException e)
			{
				throw new RuntimeException(e);
			}
		}
		return file;
	}

	public static void draw(IAtomContainer mol[], List<String> props) throws CDKException
	{
		String[] smiles = new String[mol.length];
		String[] info = new String[mol.length];
		for (int i = 0; i < smiles.length; i++)
		{
			smiles[i] = CDKConverter.toSmiles(mol[i]);
			String inf = "";
			for (String p : props)
				inf += mol[i].getProperty(p) + "<br>";
			info[i] = inf;
		}
		draw(smiles, info);
	}

	public static void draw(String... smiles)
	{
		draw(smiles, null);
	}

	public static void draw(String smiles[], String info[])
	{
		JPanel p = new JPanel(new GridLayout(3, 2, 1, 1));
		for (int i = 0; i < smiles.length; i++)
		{
			StringBuffer s = new StringBuffer();
			String smi = smiles[i];
			if (smi.length() > 20)
				smi = smi.substring(0, 18) + "..";
			s.append("<html>" + smi);
			if (info != null)
				s.append("<br>" + info[i]);
			s.append("</html>");
			ImageIcon ic = new ImageIcon(drawToFileOB(smiles[i]));
			ic = ImageLoader.getResizedImage(ic, 1.0);
			JLabel l = new JLabel(s.toString(), ic, SwingConstants.CENTER);
			l.setFont(l.getFont().deriveFont((float) 10.0));
			l.setFont(l.getFont().deriveFont(Font.PLAIN));
			p.add(l);
		}
		SwingUtil.showInDialog(new JScrollPane(p), ArrayUtil.toString(smiles));
	}

}
