package org.mg.cdklib.cfp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.mg.cdklib.CDKConverter;
import org.mg.cdklib.data.DataLoader;
import org.mg.javalib.util.CountedSet;
import org.mg.javalib.util.SetUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

public class CFPMiner extends BasicCFPMiner
{
	private static final long serialVersionUID = 6L;

	protected List<String> endpoints;
	protected HashMap<String, String> trainingInchisToEndpoint;

	transient LinkedHashMap<CFPFragment, LinkedHashSet<Integer>> fragmentToCompound_unfiltered;
	transient String[] classValues;
	transient Integer activeIdx;

	public CFPMiner(List<String> endpoints) throws CDKException
	{
		this.endpoints = endpoints;
	}

	public String[] getClassValues()
	{
		if (classValues == null)
			classValues = DataLoader.getClassValues(endpoints);
		return classValues;
	}

	public int getActiveIdx()
	{
		if (activeIdx == null)
			activeIdx = DataLoader.getActiveIdx(getClassValues());
		return activeIdx;
	}

	public List<String> getEndpoints()
	{
		return endpoints;
	}

	public void mine(List<String> smiles) throws Exception
	{
		trainingInchisToEndpoint = new HashMap<>();
		int idx = 0;
		for (String smi : smiles)
		{
			String inchi = CDKConverter.toInchi(smi);
			if (trainingInchisToEndpoint.containsKey(inchi))
				throw new IllegalStateException("no duplicates allowed!");
			trainingInchisToEndpoint.put(inchi, endpoints.get(idx));
			idx++;
		}

		super.mine(smiles);
	}

	public String getTrainingActivity(String smiles) throws CDKException
	{
		return trainingInchisToEndpoint.get(CDKConverter.toInchi(smiles));
	}

	private void applyMinFreq(Set<Integer> compoundSubset, int minFreq)
	{
		// remove with min-freq 0, would not have been mined on subset
		List<CFPFragment> fragmentToDelete = new ArrayList<CFPFragment>();
		for (CFPFragment f : fragmentToCompound.keySet())
			if (SetUtil.intersectSize(fragmentToCompound.get(f), compoundSubset) == 0)
				fragmentToDelete.add(f);
		removeFragments(fragmentToDelete);

		// remove only if larger than hashfoldize
		if (fragmentToCompound.size() <= hashfoldsize)
			return;

		// collect tuples of hash code and min-freq
		List<int[]> hf = new ArrayList<int[]>();
		for (CFPFragment f : fragmentToCompound.keySet())
			hf.add(new int[] { f.getId(), SetUtil.intersectSize(fragmentToCompound.get(f), compoundSubset) });

		// sort according to min freq, in decreasing order
		Collections.sort(hf, new Comparator<int[]>()
		{
			@Override
			public int compare(int[] o1, int[] o2)
			{
				// swap params to sort in descending order
				return Integer.valueOf(o2[1]).compareTo(Integer.valueOf(o1[1]));
			}
		});

		Set<CFPFragment> toRemove = new HashSet<>();
		// apply min-freq while size > hashfoldsize
		for (int i = hashfoldsize; i < hf.size(); i++)
		{
			if (hf.get(i)[1] < absMinFreq)
				toRemove.add(new CFPFragment(hf.get(i)[0]));
		}
		removeFragments(toRemove);

		System.err.println("applied min freq filter: " + fragmentToCompound.size());

	}

	protected static <T> long[] nominalCounts(List<T> domain, List<T> values)
	{
		CountedSet<T> countedValues = CountedSet.create(values);
		long[] counts = new long[domain.size()];
		int i = 0;
		for (T v : domain)
			counts[i++] = countedValues.getCount(v);
		return counts;
	}

	public void applyClosedSetFilter(Set<Integer> compoundSubset)
	{
		if (fragmentToCompound.size() <= hashfoldsize)
			return;
		int maxNumRemove = fragmentToCompound.size() - hashfoldsize;

		LinkedHashMap<CFPFragment, Set<Integer>> fragmentToCompoundSubset = new LinkedHashMap<>();
		LinkedHashMap<Integer, List<CFPFragment>> fragmentWithSimilarSubsets = new LinkedHashMap<>();
		for (CFPFragment f : fragmentToCompound.keySet())
		{
			@SuppressWarnings("unchecked")
			Set<Integer> hCompounds = (Set<Integer>) fragmentToCompound.get(f).clone();
			hCompounds.retainAll(compoundSubset);
			fragmentToCompoundSubset.put(f, hCompounds);

			int eq = hCompounds.hashCode();
			if (!fragmentWithSimilarSubsets.containsKey(eq))
				fragmentWithSimilarSubsets.put(eq, new ArrayList<CFPFragment>());
			fragmentWithSimilarSubsets.get(eq).add(f);
		}

		if (fragmentList == null)
			getFragmentViaIdx(0);
		Set<CFPFragment> fragmentsToRemove = new HashSet<>();

		for (List<CFPFragment> fragmentSet : fragmentWithSimilarSubsets.values())
		{
			for (int i = 0; i < fragmentSet.size() - 1; i++)
			{
				CFPFragment f1 = fragmentSet.get(i);
				if (fragmentsToRemove.contains(f1))
					continue;

				for (int j = i + 1; j < fragmentSet.size(); j++)
				{
					CFPFragment f2 = fragmentSet.get(j);
					if (fragmentsToRemove.contains(f2))
						continue;

					if (fragmentToCompoundSubset.get(f1).equals(fragmentToCompoundSubset.get(f2)))
					{
						CFPFragment obsolete = getNonClosed(f1, f2, fragmentToCompoundSubset.get(f1));
						if (obsolete != null)
						{
							fragmentsToRemove.add(obsolete);
							if (fragmentsToRemove.size() >= maxNumRemove || obsolete == f1)
								break;
						}
					}
				}
				if (fragmentsToRemove.size() >= maxNumRemove)
					break;
			}
			if (fragmentsToRemove.size() >= maxNumRemove)
				break;
		}

		removeFragments(fragmentsToRemove);

		System.err.println("applied closed fragment filter: " + fragmentToCompound.size());
	}

	public CFPFragment getNonClosed(CFPFragment f1, CFPFragment f2, Set<Integer> compounds)
	{

		//		System.out.println("check if " + h1 + " is superset of " + h2 + ", compounds: " + compounds);
		try
		{
			boolean f1SupersetCandiate = true;
			boolean f2SupersetCandiate = true;
			for (Integer c : compounds)
			{
				IAtomContainer mol = CDKConverter.parseSmiles(trainingDataSmiles.get(c));
				Set<Integer> atoms1 = getAtomsMultiple(mol, f1);
				Set<Integer> atoms2 = getAtomsMultiple(mol, f2);
				//				System.out.println("mol " + c + " atoms h1: " + atoms1 + " atoms h2: " + atoms2);
				if (f1SupersetCandiate)
					if (!SetUtil.isSubSet(atoms1, atoms2))
						f1SupersetCandiate = false;
				if (f2SupersetCandiate)
					if (!SetUtil.isSubSet(atoms2, atoms1))
						f2SupersetCandiate = false;
				if (!f1SupersetCandiate && !f2SupersetCandiate)
					return null;
			}
			if (f1SupersetCandiate)
				return f2;
			else if (f2SupersetCandiate)
				return f1;
			else
				throw new IllegalStateException();
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}

	public void applyChiSquareFilter(Set<Integer> compoundSubset)
	{
		// chi square is apply to a counts array for each class
		// e.g. it compares 45 x active 41 x inactive in the compoundSubset
		// to feature-x with 31 x active, 5 x inactive for all compounds with feature-value(feature-x, compound-y)= active
		if (endpoints.size() != numCompounds)
			throw new IllegalArgumentException(endpoints.size() + " != " + numCompounds);
		if (fragmentToCompound.size() <= hashfoldsize)
			return;
		List<String> domain = new ArrayList<String>(new HashSet<String>(endpoints));
		List<String> subsetEndpoints = new ArrayList<String>();
		for (Integer c : compoundSubset)
			subsetEndpoints.add(endpoints.get(c));
		long[] all = nominalCounts(domain, subsetEndpoints);
		//		System.out.println(domain);
		//		System.out.println(ArrayUtil.toString(all));

		// create tuples of hash code and p-value
		List<double[]> hp = new ArrayList<double[]>();
		for (CFPFragment f : fragmentToCompound.keySet())
		{
			List<String> values = new ArrayList<String>();
			for (Integer c : fragmentToCompound.get(f))
				if (compoundSubset.contains(c))
					values.add(endpoints.get(c));
			//			System.out.println(values);
			//			System.out.println(ArrayUtil.toString(nominalCounts(domain, values)));

			double p = Double.MAX_VALUE;
			if (values.size() > 0)
			{
				long sel[] = nominalCounts(domain, values);
				//System.err.println(ArrayUtil.toString(all));
				//				System.err.println(values.size());
				//System.err.println(ArrayUtil.toString(sel));
				p = TestUtils.chiSquareTestDataSetsComparison(sel, all);
				//				System.err.println("p value of hash code " + h + " is " + p);
				//				if (p > pValueThreshold)
				//					hashCodeToDelete.add(h);
			}
			hp.add(new double[] { f.getId(), p });
		}

		//sort tuples according to p-values
		Collections.sort(hp, new Comparator<double[]>()
		{
			@Override
			public int compare(double[] o1, double[] o2)
			{
				return Double.valueOf(o1[1]).compareTo(Double.valueOf(o2[1]));
			}
		});

		//sort compounds higher than hash
		Set<CFPFragment> toRemove = new HashSet<>();
		for (int i = hashfoldsize; i < hp.size(); i++)
		{
			//			System.err.println("remove fragment with p-value " + hp.get(i)[1]);
			toRemove.add(new CFPFragment((int) hp.get(i)[0]));
		}
		removeFragments(toRemove);
		if (fragmentToCompound.size() != hashfoldsize)
			throw new IllegalStateException();

		System.err.println("applied chi square filter: " + fragmentToCompound.size());
	}

	private void removeFragments(Collection<CFPFragment> fragmentToDelete)
	{
		for (CFPFragment f : fragmentToDelete)
			fragmentToCompound.remove(f);
		compoundToFragment = null;
		fragmentList = null;
	}

	protected void minePairs(Set<Integer> compoundSubset)
	{
	}

	public void applyFilter()
	{
		HashSet<Integer> allCompounds = new HashSet<Integer>();
		for (int i = 0; i < numCompounds; i++)
			allCompounds.add(i);
		applyFilter(allCompounds);
	}

	@SuppressWarnings("unchecked")
	public void applyFilter(Set<Integer> filterSubset)
	{
		if (featureSelection != FeatureSelection.filt)
			return;

		// undo old filter
		if (fragmentToCompound_unfiltered != null)
		{
			fragmentToCompound = (LinkedHashMap<CFPFragment, LinkedHashSet<Integer>>) fragmentToCompound_unfiltered
					.clone();
			compoundToFragment = null;
			fragmentList = null;
		}
		else
			fragmentToCompound_unfiltered = (LinkedHashMap<CFPFragment, LinkedHashSet<Integer>>) fragmentToCompound
					.clone();

		System.err.println("apply filtering: " + fragmentToCompound.size());

		// apply new filter
		applyMinFreq(filterSubset, absMinFreq);

		minePairs(filterSubset);
		applyClosedSetFilter(filterSubset);

		applyChiSquareFilter(filterSubset);

		//		System.out.println("filtered to: " + this);
	}
}
