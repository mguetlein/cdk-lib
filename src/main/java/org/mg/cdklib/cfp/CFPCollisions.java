package org.mg.cdklib.cfp;

import java.io.File;

import org.mg.cdklib.data.CDKDataset;
import org.mg.cdklib.data.DataProvider;
import org.mg.cdklib.data.DataProvider.DataID;
import org.mg.javalib.datamining.ResultSet;
import org.mg.javalib.datamining.ResultSetIO;

public class CFPCollisions
{

	public static void print() throws Exception
	{
		ResultSet res = new ResultSet();
		//		CFPType types[] = new CFPType[] { CFPType.ecfp6, CFPType.ecfp4, CFPType.ecfp2, CFPType.ecfp0 };
		CFPType types[] = new CFPType[] { CFPType.fcfp6, CFPType.fcfp4, CFPType.fcfp2,
				CFPType.fcfp0 };
		int dCount = 0;
		for (DataID name : DataProvider.cfpDatasets())
		{
			System.out.println(dCount + ": " + name);
			CDKDataset d = DataProvider.getDataset(name);

			for (CFPType type : types)
			{
				//			if (!name.startsWith("CPDBAS") && !name.startsWith("AMES") && !name.startsWith("NCTRER"))
				//				continue;

				System.out.println(type);

				int idx = res.addResult();

				CFPMiner miner = new CFPMiner(d.getEndpoints());
				miner.type = type;
				miner.featureSelection = FeatureSelection.filt;
				miner.hashfoldsize = 1024;
				miner.mine(d.getSmiles());

				res.setResultValue(idx, "Dataset", name);
				res.setResultValue(idx, "Type", type + "");
				res.setResultValue(idx, "Compounds", miner.getNumCompounds());
				res.setResultValue(idx, "Fragments", miner.getNumFragments());

				for (int size : new int[] { 1024, 2048, 4096, 8192 })
				{
					miner = new CFPMiner(d.getEndpoints());
					miner.type = type;
					miner.featureSelection = FeatureSelection.fold;
					miner.hashfoldsize = size;
					miner.mine(d.getSmiles());
					miner.estimateCollisions(res, idx, size + " ");
				}
			}

			res.sortResults("Dataset");
			System.out.println("\n");
			System.out.println(res.toNiceString());

			if (types.length > 1)
			{
				ResultSet joined = res.copy().join("Type");
				joined.removePropery("Dataset");

				System.out.println("\n");
				System.out.println(joined.toNiceString());
			}

			dCount++;
			//			if (dCount > 2)
			//				break;
		}
		ResultSetIO
				.printToTxtFile(
						new File(System.getProperty("user.home")
								+ "/results/cdklib/data_collisions/collisions_fcfp.result"),
						res, true);
		System.exit(1);
	}

	public static void main(String[] args) throws Exception
	{
		print();
	}
}
