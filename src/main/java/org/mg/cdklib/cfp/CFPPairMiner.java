package org.mg.cdklib.cfp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.mg.cdklib.CDKConverter;
import org.mg.cdklib.data.CDKDataset;
import org.mg.cdklib.data.DataProvider;
import org.mg.cdklib.data.DataProvider.Dataset;
import org.mg.javalib.util.CountedSet;
import org.mg.javalib.util.HashUtil;
import org.mg.javalib.util.ListUtil;
import org.mg.javalib.util.StringUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

public class CFPPairMiner extends CFPMiner
{
	private static final long serialVersionUID = 1L;

	private transient List<Pair> pairs;

	public CFPPairMiner(List<String> endpoints) throws CDKException
	{
		super(endpoints);
	}

	public class Pair
	{
		CFPFragment f1;
		CFPFragment f2;
		Set<Integer> commonCompounds;
		Set<Integer> commonAdjacentCompounds;
		double pCommon;
		double pCommonAdj;
		double pDiff;

		public Pair(CFPFragment f1, CFPFragment f2, Set<Integer> commonCompounds,
				Set<Integer> commonAdjacentCompounds)
		{
			if (f2.getId() < f1.getId())
			{
				this.f1 = f2;
				this.f2 = f1;
			}
			else
			{
				this.f1 = f1;
				this.f2 = f2;
			}
			this.commonCompounds = commonCompounds;
			this.commonAdjacentCompounds = commonAdjacentCompounds;

			pCommonAdj = pValue(commonAdjacentCompounds);
			pCommon = pValue(commonCompounds);
			pDiff = pCommon - pCommonAdj;
		}

		public boolean equals(Object o)
		{
			return (o instanceof Pair && f1 == ((Pair) o).f1 && f2 == ((Pair) o).f2);
		}

		@Override
		public int hashCode()
		{
			return HashUtil.hashCode(f1, f2);
		}

		@Override
		public String toString()
		{
			String str = "pDiff:" + StringUtil.formatDouble(pDiff) + " ";
			List<String> values = new ArrayList<String>();
			for (Integer c : commonCompounds)
				values.add(endpoints.get(c));
			str += "common:" + CountedSet.create(values) + "(p:" + StringUtil.formatDouble(pCommon)
					+ ") ";
			values = new ArrayList<String>();
			for (Integer c : commonAdjacentCompounds)
				values.add(endpoints.get(c));
			str += "adjacent:" + CountedSet.create(values) + "(p:"
					+ StringUtil.formatDouble(pCommonAdj) + ") ";
			return str;
		}

		public String getName()
		{
			return "pair-" + f1 + "-" + f2;
		}

		public boolean isCompoundAdjacent(Integer compound)
		{
			return commonAdjacentCompounds.contains(compound);
		}

	}

	transient List<String> domain;
	transient long[] all;

	private double pValue(Set<Integer> compounds)
	{
		List<String> values = new ArrayList<String>();
		for (Integer c : compounds)
			values.add(endpoints.get(c));
		long sel[] = nominalCounts(domain, values);
		return TestUtils.chiSquareTestDataSetsComparison(sel, all);
	}

	@Override
	public String getName()
	{
		return super.getName() + "_pairs";
	}

	@Override
	@SuppressWarnings("unchecked")
	protected void minePairs(Set<Integer> compoundSubset)
	{
		if (fragmentList == null)
			getFragmentViaIdx(0);

		domain = new ArrayList<String>(new HashSet<String>(endpoints));
		List<String> subsetEndpoints = new ArrayList<String>();
		for (Integer c : compoundSubset)
		{
			//			System.out.println(c);
			subsetEndpoints.add(endpoints.get(c));
		}
		System.out.println("subset endpoints " + CountedSet.create(subsetEndpoints));
		all = nominalCounts(domain, subsetEndpoints);

		pairs = new ArrayList<>();

		int minFrequency = Math.max(1, (int) (compoundSubset.size() * 0.05));
		System.out.println("min-freq for pairs: " + minFrequency);
		for (int i = 0; i < fragmentList.length - 1; i++)
		{
			if (i % 100 == 0)
				System.out.println(i);

			CFPFragment f1 = fragmentList[i];
			if (fragmentToCompound.get(f1).size() < minFrequency)
				continue;
			Set<Integer> h1Compounds = (Set<Integer>) fragmentToCompound.get(f1).clone();
			h1Compounds.retainAll(compoundSubset);
			if (h1Compounds.size() < minFrequency)
				continue;

			for (int j = i + 1; j < fragmentList.length; j++)
			{
				CFPFragment f2 = fragmentList[j];
				if (fragmentToCompound.get(f2).size() < minFrequency)
					continue;
				HashSet<Integer> intersect = new HashSet<Integer>(h1Compounds);
				intersect.retainAll(fragmentToCompound.get(f2));
				if (intersect.size() < minFrequency)
					continue;
				Set<Integer> adj = getAdjacent(f1, f2, intersect);
				if (adj == null) // overlap
					continue;
				if (adj.size() < minFrequency)
					continue;
				//				System.out.println(numAdj + "/" + intersect.size());
				Pair pair = new Pair(f1, f2, intersect, adj);
				//				pairAdjacent.put(new Pair(h1, h2), adj);
				if (pair.pDiff > 0.1)
					pairs.add(pair);
			}
		}

		Collections.sort(pairs, new Comparator<Pair>()
		{
			@Override
			public int compare(Pair o1, Pair o2)
			{
				return Double.valueOf(o2.pDiff).compareTo(o1.pDiff);
			}
		});
		int idx = 0;
		for (Pair pair : pairs)
		{
			System.out.println((idx++) + " " + pair);
		}
	}

	private Set<Integer> getAdjacent(CFPFragment h1, CFPFragment h2,
			HashSet<Integer> intersectCompounds)
	{
		try
		{
			Set<Integer> adj = new HashSet<>();
			for (Integer c : intersectCompounds)
			{
				IAtomContainer mol = CDKConverter.parseSmiles(trainingDataSmiles.get(c));
				//			int[] atoms1 = getAtoms(mol, h1);
				//			int[] atoms2 = getAtoms(mol, h2);
				Set<Integer> atoms1 = getAtomsMultiple(mol, h1);
				Set<Integer> atoms2 = getAtomsMultiple(mol, h2);

				boolean adjacent = false;
				for (int a1 : atoms1)
				{
					List<IAtom> connected = null;
					if (!adjacent)
						connected = mol.getConnectedAtomsList(mol.getAtom(a1));

					for (int a2 : atoms2)
					{
						if (a1 == a2)
							return null;
						if (!adjacent && connected.contains(mol.getAtom(a2)))
							adjacent = true;
					}
				}
				if (adjacent)
					adj.add(c);
			}
			return adj;
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}

	public static void demo() throws Exception
	{
		Dataset d = Dataset.NCTRER;
		CFPType type = CFPType.ecfp2;
		FeatureSelection featureSelection = FeatureSelection.filt;
		int hashfoldsize = 1024;

		CDKDataset dataset = DataProvider.getDataset(d);
		List<String> list = dataset.getSmiles();
		List<String> endpointValues = dataset.getEndpoints();
		ListUtil.scramble(new Random(1), list, endpointValues);

		CFPMiner cfps = new CFPPairMiner(ListUtil.cast(String.class, endpointValues));
		cfps.type = type;
		cfps.featureSelection = featureSelection;
		cfps.hashfoldsize = hashfoldsize;
		cfps.mine(list);
		System.out.println(cfps);

		cfps.applyFilter();
	}

	public static void main(String[] args) throws Exception
	{
		demo();
	}

	public int getNumPairs()
	{
		return pairs.size();
	}

	public Pair getPair(int idx)
	{
		return pairs.get(idx);
	}
}
