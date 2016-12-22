package org.mg.cdklib.match;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.BitSet;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.SerializationUtils;
import org.openscience.cdk.exception.CDKException;

public class PersistentSubstructureCountMatcher implements SubstructureCountMatcher
{
	DefaultSubstructureCountMatcher scm;
	String file;
	int lastStoreCount;
	private static int STORE_INTERVAL = 100;

	PersistentSubstructureCountMatcher(String smarts[], String file)
	{
		this.file = file;
		if (new File(file).exists())
		{
			try
			{
				scm = SerializationUtils.deserialize(new FileInputStream(file));
				lastStoreCount = scm.numMatchedSmiles();
				System.err.println(
						"substructure-matcher> read '" + scm.numMatchedSmiles() + "' from " + file);
			}
			catch (FileNotFoundException e)
			{
				throw new RuntimeException(e);
			}
		}
		else
		{
			scm = new DefaultSubstructureCountMatcher(smarts);
		}
	}

	public void store()
	{
		if (scm.numMatchedSmiles() == lastStoreCount)
			return;
		try
		{
			System.err.println(
					"substructure-matcher> storing '" + scm.numMatchedSmiles() + "' to " + file);
			File tmp = File.createTempFile("count", "matcher");
			SerializationUtils.serialize(scm, new FileOutputStream(tmp));
			FileUtils.deleteQuietly(new File(file));
			FileUtils.moveFile(tmp, new File(file));
			lastStoreCount = scm.numMatchedSmiles();
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

	public short[] matchCount(String smiles) throws CDKException
	{
		try
		{
			return scm.matchCount(smiles);
		}
		finally
		{
			if (scm.numMatchedSmiles() % STORE_INTERVAL == 0)
				store();
		}
	}

	@Override
	public BitSet match(String smiles) throws CDKException
	{
		try
		{
			return scm.match(smiles);
		}
		finally
		{
			if (scm.numMatchedSmiles() % STORE_INTERVAL == 0)
				store();
		}
	}

	@Override
	public String getSmarts(int idx)
	{
		return scm.getSmarts(idx);
	}

	@Override
	public int size()
	{
		return scm.size();
	}
}
