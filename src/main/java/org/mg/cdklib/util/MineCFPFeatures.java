package org.mg.cdklib.util;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.mg.cdklib.cfp.BasicCFPMiner;
import org.mg.cdklib.cfp.CFPFragment;
import org.mg.cdklib.cfp.CFPType;
import org.mg.cdklib.cfp.FeatureSelection;
import org.mg.cdklib.data.CDKDataset;
import org.mg.cdklib.data.DataLoader;
import org.mg.javalib.util.CollectionUtil;

public class MineCFPFeatures
{
	public static void main(String[] args) throws Exception
	{

		CDKDataset d = DataLoader.getDatasetFromSMILES("mastermatrix",
				"/home/martin/data/envipath/mastermatrix.ob.smi");

		BasicCFPMiner miner = new BasicCFPMiner();
		miner.setType(CFPType.ecfp4);
		miner.setFeatureSelection(FeatureSelection.none);
		miner.mine(d.getSmiles());
		System.out.println(miner);

		List<String> header = new ArrayList<>();
		header.add("Cpid");
		header.add("Smiles");
		for (int i = 0; i < miner.getNumFragments(); i++)
			header.add(miner.getFragmentViaIdx(i).getId() + "");

		CSVPrinter csv = new CSVPrinter(
				new FileWriter(new File("/home/martin/data/envipath/mastermatrix.features."
						+ miner.getCFPType() + ".csv")),
				CSVFormat.DEFAULT.withHeader(CollectionUtil.toArray(header)));

		for (int i = 0; i < miner.getNumCompounds(); i++)
		{
			List<Object> vals = new ArrayList<>();
			vals.add(d.getEndpoints().get(i));
			if (!d.getSmiles().get(i).equals(miner.getTrainingDataSmiles().get(i)))
				throw new IllegalArgumentException();
			vals.add(d.getSmiles().get(i));
			Set<CFPFragment> frags = miner.getFragmentsForCompound(i);
			for (int j = 0; j < miner.getNumFragments(); j++)
				vals.add(frags.contains(miner.getFragmentViaIdx(j)) ? 1 : 0);
			csv.printRecord(vals);
		}

		csv.close();

	}
}
