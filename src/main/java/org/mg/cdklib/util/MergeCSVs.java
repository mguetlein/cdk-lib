package org.mg.cdklib.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.mg.cdklib.CDKConverter;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.DoubleArraySummary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

public class MergeCSVs
{
	static HashMap<String, HashMap<String, Object>> compounds = new HashMap<>();

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException, InvalidSmilesException, CDKException
	{
		String n[] = { "swiss", "mazzatorta" };
		for (String s : n)
		{
			CSVParser csv = new CSVParser(
					new FileReader(new File("/home/martin/data/nestle6/" + s + ".csv")),
					CSVFormat.DEFAULT.withHeader());

			System.out.println(csv.getHeaderMap());
			int i = 0;
			for (CSVRecord row : csv)
			{
				String origSmiles = row.get("SMILES");
				String adjustedSmiles;
				if (origSmiles.equals("Cn1n(C)c(cc1c1ccccc1)c1ccccc1"))
					adjustedSmiles = "CN1N(C)C(CC1c1ccccc1)c1ccccc1";
				else
					adjustedSmiles = null;

				String absSmiles = CDKConverter
						.toAbsoluteSmiles(adjustedSmiles != null ? adjustedSmiles : origSmiles);
				Double val = Double.valueOf(row.get("LOAEL"));
				String dataset = row.get("Dataset");
				if (!dataset.equals(s))
					throw new IllegalStateException();
				int idx = i++;
				//				System.out.println(dataset);

				if (!compounds.containsKey(absSmiles))
				{
					compounds.put(absSmiles, new HashMap<String, Object>());
					compounds.get(absSmiles).put("LOAEL", new ArrayList<Double>());
				}

				compounds.get(absSmiles).put(dataset, true);
				((List<Double>) compounds.get(absSmiles).get("LOAEL")).add(val);
			}

			csv.close();
			System.out.println(compounds.size());
		}

		for (String absSmiles : compounds.keySet())
		{
			String d = null;
			for (String s : n)
				if (compounds.get(absSmiles).containsKey(s))
					if (d == null)
						d = s;
					else
						d = "both";
			compounds.get(absSmiles).put("dataset", d);
		}

		List<String[]> out = new ArrayList<>();
		String header[] = { "SMILES", "dataset", "swiss", "mazzatorta", "LOAEL", "Log-LOAEL",
				"measurements" };
		out.add(header);
		for (String absSmiles : compounds.keySet())
		{
			//			System.out.println(absSmiles);
			HashMap<String, Object> vals = compounds.get(absSmiles);
			List<Object> row = new ArrayList<>();
			row.add(absSmiles);
			row.add(vals.get("dataset"));
			row.add(vals.get("swiss"));
			row.add(vals.get("mazzatorta"));
			List<Double> v = (List<Double>) vals.get("LOAEL");
			row.add(DoubleArraySummary.create(v).getMean());
			row.add(Math.log(DoubleArraySummary.create(v).getMean()));
			row.add(v.size());
			out.add(ArrayUtil.toStringArray(ArrayUtil.toArray(Object.class, row)));
		}

		BufferedWriter w = new BufferedWriter(new FileWriter("/home/martin/data/nestle6/uniq.csv"));
		for (String row[] : out)
		{
			w.write(ArrayUtil.toCSVString(row));
			w.write('\n');
		}
		w.close();
	}
}
