package org.mg.cdklib.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.io.Charsets;
import org.apache.commons.io.IOUtils;
import org.codehaus.jettison.json.JSONArray;
import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;

public class PubChemSampler
{
	int maxCID = 119000000;

	Random random = new Random();

	boolean debug = true;

	public int randomCID()
	{
		return random.nextInt(maxCID);
	}

	public String smilesFromCid(int cid)
	{
		try
		{
			String url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/" + cid + "/JSON";
			JSONObject obj = new JSONObject(IOUtils.toString(new URL(url), Charsets.UTF_8));
			obj = obj.getJSONArray("PC_Compounds").getJSONObject(0);
			JSONArray props = obj.getJSONArray("props");
			for (int i = 0; i < props.length(); i++)
				if (props.getJSONObject(i).getJSONObject("urn").getString("label").equals("SMILES")
						&& props.getJSONObject(i).getJSONObject("urn").getString("name")
								.equals("Canonical"))
					return props.getJSONObject(i).getJSONObject("value").getString("sval");
			return null;
		}
		catch (Exception e)
		{
			System.err.println("smiles for cid '" + cid + "' failed");
			return null;
		}
	}

	public HashMap<Integer, String> sample(int n)
	{
		HashMap<Integer, String> map = new HashMap<>();
		for (int i = 0; i < n; i++)
		{
			String smiles = null;
			int cid = randomCID();
			while (smiles == null)
			{
				while (map.containsKey(cid))
					cid = randomCID();
				smiles = smilesFromCid(cid);
			}
			if (debug)
				System.out
						.println("[" + i + "] sampled cid: '" + cid + "' smiles: '" + smiles + "'");
			map.put(cid, smiles);
		}
		return map;
	}

	public void sampleToSmilesFile(String file, int n)
	{
		BufferedWriter buffy = null;
		try
		{
			buffy = new BufferedWriter(new FileWriter(file));
			Set<Integer> cids = new HashSet<>();
			for (int i = 0; i < n; i++)
			{
				String smiles = null;
				int cid = randomCID();
				while (smiles == null)
				{
					while (cids.contains(cid))
						cid = randomCID();
					smiles = smilesFromCid(cid);
				}
				if (debug)
					System.out.println(
							"[" + i + "] sampled cid: '" + cid + "' smiles: '" + smiles + "'");
				buffy.write(smiles + "\t" + cid + "\n");
				if (i % 100 == 0)
					buffy.flush();
			}
			if (debug)
				System.out.println("stored '" + n + "' sampled compounds to '" + file + "'");
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
		finally
		{
			IOUtils.closeQuietly(buffy);
		}
	}

	public static void main(String[] args) throws MalformedURLException, JSONException, IOException
	{
		PubChemSampler s = new PubChemSampler();
		//		HashMap<Integer, String> map = s.sample(3);
		//		for (Integer cid : map.keySet())
		//		{
		//			System.out.println(cid + " " + map.get(cid));
		//		}
		s.sampleToSmilesFile("/tmp/pubchem.smi", 3000);
	}
}
