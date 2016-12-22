package org.mg.cdklib.cfp;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.mg.cdklib.CDKConverter;
import org.mg.cdklib.data.CDKDataset;
import org.mg.cdklib.data.DataLoader;
import org.mg.javalib.datamining.ResultSet;
import org.mg.javalib.util.CountedSet;
import org.mg.javalib.util.DoubleArraySummary;
import org.mg.javalib.util.FileUtil;
import org.mg.javalib.util.HashUtil;
import org.mg.javalib.util.SetUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;

public class BasicCFPMiner implements Serializable
{
	private static final long serialVersionUID = 7L;

	protected int numCompounds = 0;
	protected List<String> trainingDataSmiles;
	protected LinkedHashMap<CFPFragment, LinkedHashSet<Integer>> fragmentToCompound = new LinkedHashMap<CFPFragment, LinkedHashSet<Integer>>();
	protected int numUnfoldedConflicts = 0;

	protected CFPType type;
	protected FeatureSelection featureSelection;
	protected int hashfoldsize;
	protected int absMinFreq = 2;

	transient CircularFingerprinter fp;
	transient CFPFragment[] fragmentList;
	transient HashMap<Integer, LinkedHashSet<CFPFragment>> compoundToFragment;
	transient HashMap<CFPFragment, Integer> fragmentToIdx = new HashMap<CFPFragment, Integer>();
	transient HashMap<CFPFragment, Integer> fragmentToIteration = new HashMap<CFPFragment, Integer>();
	transient HashMap<CFPFragment, Integer> fragmentToNumAtoms = new HashMap<CFPFragment, Integer>();
	transient HashMap<IAtomContainer, LinkedHashSet<CFPFragment>> testMoleculeToFragment;

	transient HashMap<CFPFragment, LinkedHashSet<CFPFragment>> subFragments;
	transient HashMap<CFPFragment, LinkedHashSet<CFPFragment>> superFragments;

	private HashMap<Integer, Set<Integer>> collisionMap;

	public BasicCFPMiner()
	{
	}

	public List<String> getTrainingDataSmiles()
	{
		return trainingDataSmiles;
	}

	public int getHashfoldsize()
	{
		return hashfoldsize;
	}

	public int getNumFragments()
	{
		return fragmentToCompound.size();
	}

	public CFPType getCFPType()
	{
		return type;
	}

	public FeatureSelection getFeatureSelection()
	{
		return featureSelection;
	}

	public String getNiceFragmentDescription()
	{
		return featureSelection.attribute() + " " + type.toNiceString() + " fragments";
	}

	public void setType(CFPType type)
	{
		this.type = type;
	}

	public void setFeatureSelection(FeatureSelection featureSelection)
	{
		this.featureSelection = featureSelection;
	}

	public void setHashfoldsize(int hashfoldsize)
	{
		this.hashfoldsize = hashfoldsize;
	}

	public void setAbsMinFreq(int absMinFreq)
	{
		this.absMinFreq = absMinFreq;
	}

	public int getAbsMinFreq()
	{
		return absMinFreq;
	}

	public boolean isFragmentIncludedInCompound(int compound, CFPFragment fragment)
	{
		if (fragmentToCompound.get(fragment) == null)
			throw new IllegalStateException("no compounds for fragment, should have been removed! "
					+ fragment + " " + fragmentToCompound.get(fragment));
		return (fragmentToCompound.get(fragment).contains(compound));
	}

	private void initCircularFingerprinter()
	{
		collisionMap = new HashMap<>();
		fp = new CircularFingerprinter(type.getClassType())
		{
			public IBitFingerprint getBitFingerprint(IAtomContainer mol) throws CDKException
			{
				calculate(mol);
				if (featureSelection != FeatureSelection.fold)
					return null;
				final BitSet bits = new BitSet(hashfoldsize);
				for (int n = 0; n < getFPCount(); n++)
				{
					int i = getFP(n).hashCode;
					long b = i >= 0 ? i : ((i & 0x7FFFFFFF) | (1L << 31));
					int bit = (int) (b % hashfoldsize);
					if (!collisionMap.containsKey(bit))
						collisionMap.put(bit, new HashSet<Integer>());
					collisionMap.get(bit).add(i);
					bits.set(bit);
				}
				return new BitSetFingerprint(bits);
			}
		};
	}

	public int[] getAtoms(String smiles, CFPFragment fragment) throws CDKException
	{
		return getAtoms(CDKConverter.parseSmiles(smiles), fragment);
	}

	public int[] getAtoms(IAtomContainer mol, CFPFragment fragment) throws CDKException
	{
		if (featureSelection == FeatureSelection.fold)
			throw new IllegalArgumentException();
		int atoms[] = null;
		initCircularFingerprinter();
		fp.getBitFingerprint(mol);
		for (int i = 0; i < fp.getFPCount(); i++)
		{
			if (fp.getFP(i).hashCode == fragment.getId())
			{
				atoms = fp.getFP(i).atoms;
				break;
			}
		}
		return atoms;
	}

	private static transient HashMap<Integer, Set<Integer>> atomsMultCache = new HashMap<>();

	/**
	 * fragment may occur multiple times
	 * 
	 * @param mol
	 * @param fragment
	 * @return
	 * @throws Exception
	 */
	public Set<Integer> getAtomsMultiple(IAtomContainer mol, CFPFragment fragment)
			throws CDKException
	{
		if (featureSelection == FeatureSelection.fold)
			throw new IllegalArgumentException();

		Integer key = HashUtil.hashCode(type, mol, fragment);
		if (!atomsMultCache.containsKey(key))
		{
			Set<Integer> atoms = new HashSet<>();
			initCircularFingerprinter();
			fp.getBitFingerprint(mol);
			for (int i = 0; i < fp.getFPCount(); i++)
				if (fp.getFP(i).hashCode == fragment.getId())
					for (int a : fp.getFP(i).atoms)
						atoms.add(a);
			atomsMultCache.put(key, atoms);
		}
		return atomsMultCache.get(key);
	}

	public static void main(String[] args) throws Exception
	{
		CDKDataset d = DataLoader.INSTANCE.getDataset("AMES");
		CFPMiner miner = new CFPMiner(d.getEndpoints());
		miner.type = CFPType.ecfp4;
		miner.featureSelection = FeatureSelection.filt;
		miner.hashfoldsize = 4096;
		miner.mine(d.getSmiles());
		miner.applyFilter();

		System.out.println(miner.toString());

		for (Integer fragIdx : new Integer[] { 4, 311 })
		{
			CFPFragment f = miner.getFragmentViaIdx(fragIdx);
			System.out.println(fragIdx + ": " + f.getId());
			System.out.println("matches " + miner.getCompoundsForFragment(f).size());
		}
		//		miner.mineSubAndSuperFragments();

	}

	private static transient HashMap<Integer, Set<Set<Integer>>> atomsMultCacheDistinct = new HashMap<>();

	public Set<Set<Integer>> getAtomsMultipleDistinct(IAtomContainer mol, CFPFragment fragment)
			throws CDKException
	{
		if (featureSelection == FeatureSelection.fold)
			throw new IllegalArgumentException();

		Integer key = HashUtil.hashCode(type, mol, fragment);
		if (!atomsMultCacheDistinct.containsKey(key))
		{
			Set<Set<Integer>> atomsDistinct = new HashSet<>();
			initCircularFingerprinter();
			fp.getBitFingerprint(mol);
			for (int i = 0; i < fp.getFPCount(); i++)
				if (fp.getFP(i).hashCode == fragment.getId())
				{
					Set<Integer> atoms = new HashSet<>();
					for (int a : fp.getFP(i).atoms)
						atoms.add(a);
					atomsDistinct.add(atoms);
				}
			atomsMultCacheDistinct.put(key, atomsDistinct);
		}
		return atomsMultCacheDistinct.get(key);
	}

	public void mine(List<String> smiles) throws Exception
	{
		this.trainingDataSmiles = smiles;

		initCircularFingerprinter();
		for (String smi : smiles)
		{
			IBitFingerprint finger = fp.getBitFingerprint(CDKConverter.parseSmiles(smi));
			if (featureSelection == FeatureSelection.fold)
				for (int i : finger.getSetbits())
					insert(fragmentToCompound, new CFPFragment(i), numCompounds);
			else
				for (int i = 0; i < fp.getFPCount(); i++)
				{
					//					if (fp.getFP(i).hashCode == 1039376976)
					//					{
					//						System.err.println(fp.getFP(i).atoms.length + " " + ArrayUtil.toString(fp.getFP(i).atoms));
					//						try
					//						{
					//							drawFP(null, mol, fp.getFP(i).atoms);
					//						}
					//						catch (Exception e)
					//						{
					//							e.printStackTrace();
					//						}
					//					}
					CFPFragment frag = new CFPFragment(fp.getFP(i).hashCode);
					insert(fragmentToCompound, frag, numCompounds);
					boolean conflict = check(fragmentToIteration, frag, fp.getFP(i).iteration);
					conflict |= check(fragmentToNumAtoms, frag, fp.getFP(i).atoms.length);
					if (conflict)
					{
						//						System.err.println(fp.getFP(i).hashCode);
						//						if (numUnfoldedConflicts > 5)
						//							System.exit(1);
						numUnfoldedConflicts++;
					}
				}
			numCompounds++;
		}
	}

	private static <T1, T2> void insert(HashMap<T1, LinkedHashSet<T2>> map, T1 key, T2 val)
	{
		if (!map.containsKey(key))
			map.put(key, new LinkedHashSet<T2>());
		map.get(key).add(val);
	}

	private static <T> boolean check(HashMap<T, Integer> map, T key, int val)
	{
		if (map.containsKey(key) && map.get(key) != val)
			return true; //System.err.println("conflict " + key + " val1: " + val + " val2:" + map.get(key));
		map.put(key, val);
		return false;
	}

	public LinkedHashSet<Integer> getCompoundsForFragment(CFPFragment fragment)
	{
		return fragmentToCompound.get(fragment);
	}

	public LinkedHashSet<CFPFragment> getFragmentsForCompound(Integer compound)
	{
		if (compoundToFragment == null)
		{
			compoundToFragment = new HashMap<Integer, LinkedHashSet<CFPFragment>>();
			for (CFPFragment f : fragmentToCompound.keySet())
				for (Integer c : fragmentToCompound.get(f))
					insert(compoundToFragment, c, f);
		}
		if (compoundToFragment.containsKey(compound))
			return compoundToFragment.get(compound);
		else
			return new LinkedHashSet<CFPFragment>();
	}

	public double getTanimotoSimilarity(int i, int j)
	{
		HashSet<CFPFragment> h1 = getFragmentsForCompound(i);
		HashSet<CFPFragment> h2 = getFragmentsForCompound(j);
		int and = SetUtil.intersectSize(h1, h2);
		int or = h1.size() + h2.size() - and;
		return and / (double) or;
	}

	public CFPFragment getFragmentViaIdx(int fragmentIdx)
	{
		if (fragmentList == null)
		{
			fragmentList = new CFPFragment[fragmentToCompound.size()];
			int idx = 0;
			for (CFPFragment h : fragmentToCompound.keySet())
				fragmentList[idx++] = h;
		}
		return fragmentList[fragmentIdx];
	}

	public Integer getIdxForFragment(CFPFragment frag)
	{
		if (fragmentToIdx == null)
		{
			fragmentToIdx = new HashMap<>();
			int idx = 0;
			for (CFPFragment h : fragmentToCompound.keySet())
				fragmentToIdx.put(h, idx++);
		}
		return fragmentToIdx.get(frag);
	}

	transient HashMap<String, HashMap<CFPFragment, LinkedHashSet<CFPFragment>>> includedFragments;

	public Set<CFPFragment> getIncludedFragments(CFPFragment f, String smiles) throws CDKException
	{
		if (includedFragments == null || !includedFragments.containsKey(smiles))
			mineIncludedFragments(smiles);
		return includedFragments.get(smiles).get(f);
	}

	private void mineIncludedFragments(String smiles) throws CDKException
	{
		HashMap<CFPFragment, LinkedHashSet<CFPFragment>> map = new HashMap<>();

		IAtomContainer mol = CDKConverter.parseSmiles(smiles);
		List<CFPFragment> frags = new ArrayList<>();
		for (CFPFragment f : getFragmentsForTestCompound(mol))
			if (fragmentToCompound.containsKey(f))
				frags.add(f);

		for (int i1 = 0; i1 < frags.size() - 1; i1++)
		{
			CFPFragment f1 = frags.get(i1);
			Set<Integer> a1 = getAtomsMultiple(mol, f1);

			for (int i2 = i1 + 1; i2 < frags.size(); i2++)
			{
				CFPFragment f2 = frags.get(i2);
				Set<Integer> a2 = getAtomsMultiple(mol, f2);

				if (SetUtil.isSubSet(a1, a2))
					insert(map, f1, f2);
				else if (SetUtil.isSubSet(a2, a1))
					insert(map, f2, f1);
			}
		}
		//		for (int i = 0; i < frags.size(); i++)
		//			System.out.println("included in " + frags.get(i) + " : " + map.get(frags.get(i)));

		if (includedFragments == null)
			includedFragments = new HashMap<>();
		includedFragments.put(smiles, map);
	}

	public Set<CFPFragment> getSubFragments(CFPFragment frag) throws CDKException
	{
		if (subFragments == null)
			mineSubAndSuperFragments();
		return subFragments.get(frag);
	}

	public Set<CFPFragment> getSuperFragments(CFPFragment frag) throws CDKException
	{
		if (superFragments == null)
			mineSubAndSuperFragments();
		return superFragments.get(frag);
	}

	private void mineSubAndSuperFragments() throws CDKException
	{
		if (subFragments != null)
			throw new IllegalArgumentException();
		subFragments = new HashMap<CFPFragment, LinkedHashSet<CFPFragment>>();
		superFragments = new HashMap<CFPFragment, LinkedHashSet<CFPFragment>>();

		for (int i = 0; i < getNumFragments(); i++)
		{
			// get atom-matches of fragment in arbitrary mol where this fragment matches
			// (we are looking for sub-fragments of this fragment, they are in each mol) 
			CFPFragment f1 = getFragmentViaIdx(i);
			Integer mol = getCompoundsForFragment(f1).iterator().next();
			IAtomContainer molC = CDKConverter.parseSmiles(trainingDataSmiles.get(mol));

			//Set<Integer> atoms1 = getAtomsMultiple(molC, f1);
			Set<Set<Integer>> atoms1 = getAtomsMultipleDistinct(molC, f1);
			int numAtoms1 = atoms1.iterator().next().size();
			int numMatches1 = atoms1.size();

			// get all other fragments that match this mol
			for (CFPFragment f2 : getFragmentsForCompound(mol))
			{
				if (f1.equals(f2))
					continue;
				// get atom-matches for other fragment
				//Set<Integer> atoms2 = getAtomsMultiple(molC, f2);
				Set<Set<Integer>> atoms2 = getAtomsMultipleDistinct(molC, f2);
				int numAtoms2 = atoms2.iterator().next().size();
				int numMatches2 = atoms2.size();

				if (numAtoms2 < numAtoms1 && numMatches2 >= numMatches1)
				{
					boolean included2in1 = true;
					for (Set<Integer> subSet1 : atoms1)
					{
						boolean found = false;
						for (Set<Integer> subSet2 : atoms2)
						{
							if (SetUtil.isSubSet(subSet1, subSet2))
							{
								found = true;
								break;
							}
						}
						if (!found)
						{
							included2in1 = false;
							break;
						}
					}
					if (included2in1)
					{
						insert(subFragments, f1, f2);
						insert(superFragments, f2, f1);
					}
				}
			}
		}
		//		System.out.println("mined sub and super");
		//		for (int i = 0; i < getNumFragments() - 1; i++)
		//		{
		//			System.out.println("\n\n" + (i + 1));
		//			System.out.println("sub fragments:");
		//			if (subFragments.get(getFragmentViaIdx(i)) != null)
		//				for (CFPFragment f : subFragments.get(getFragmentViaIdx(i)))
		//					System.out.print((getIdxForFragment(f) + 1) + " ");
		//			System.out.println("\nsuper fragments:");
		//			if (superFragments.get(getFragmentViaIdx(i)) != null)
		//				for (CFPFragment f : superFragments.get(getFragmentViaIdx(i)))
		//					System.out.print((getIdxForFragment(f) + 1) + " ");
		//		}
		//		System.out.println("\n");
	}

	public LinkedHashSet<CFPFragment> getFragmentsForTestCompound(String smiles) throws CDKException
	{
		return getFragmentsForTestCompound(CDKConverter.parseSmiles(smiles));
	}

	public LinkedHashSet<CFPFragment> getFragmentsForTestCompound(IAtomContainer testMol)
			throws CDKException
	{
		if (testMoleculeToFragment == null)
			testMoleculeToFragment = new HashMap<IAtomContainer, LinkedHashSet<CFPFragment>>();
		if (!testMoleculeToFragment.containsKey(testMol))
		{
			if (fp == null)
				initCircularFingerprinter();
			LinkedHashSet<CFPFragment> fragments = new LinkedHashSet<CFPFragment>();
			IBitFingerprint finger = fp.getBitFingerprint(testMol);
			if (featureSelection == FeatureSelection.fold)
				for (int i : finger.getSetbits())
					fragments.add(new CFPFragment(i));
			else
				for (int i = 0; i < fp.getFPCount(); i++)
				{
					CFPFragment frag = new CFPFragment(fp.getFP(i).hashCode);
					if (fragmentToCompound.containsKey(frag))
						fragments.add(frag);
				}
			testMoleculeToFragment.put(testMol, fragments);
		}
		return testMoleculeToFragment.get(testMol);
	}

	@SuppressWarnings("unchecked")
	public BasicCFPMiner clone()
	{
		BasicCFPMiner f = new BasicCFPMiner();
		f.type = type;
		f.hashfoldsize = hashfoldsize;
		f.featureSelection = featureSelection;
		f.absMinFreq = absMinFreq;
		f.numCompounds = numCompounds;
		for (CFPFragment frag : fragmentToCompound.keySet())
			f.fragmentToCompound.put(frag,
					(LinkedHashSet<Integer>) fragmentToCompound.get(frag).clone());
		return f;
	}

	public int getNumCompounds()
	{
		return numCompounds;
	}

	public ResultSet getSummary(boolean nice)
	{
		ResultSet set = new ResultSet();
		int idx = set.addResult();
		//set.setResultValue(idx, "name", getName());
		//		if (!nice)
		//			set.setResultValue(idx, "Endpoints", CountedSet.create(endpoints));
		set.setResultValue(idx, "Num fragments", fragmentToCompound.size());
		if (!nice)
			set.setResultValue(idx, "Num compounds", numCompounds);
		set.setResultValue(idx, "Fragment type", nice ? type.toNiceString() : type);
		set.setResultValue(idx, "Feature selection",
				nice ? featureSelection.toNiceString() : featureSelection);
		//		if (featureSelection == FeatureSelection.filter)
		//		{
		//			set.setResultValue(idx, "Min frequency (relative/absolute)",
		//					relMinFreq + " / " + (int) Math.round(relMinFreq * numCompounds));
		//			set.setResultValue(idx, Character.toString((char) 967) + Character.toString((char) 178) + " max p-value",
		//					pValueThreshold);
		//		}
		if (!nice)
			set.setResultValue(idx, "Fingerprint size", hashfoldsize);
		if (!nice)
		{
			set.setResultValue(idx, "num unfolded conflicts", numUnfoldedConflicts);
			int n = 0;
			for (int c = 0; c < numCompounds; c++)
				if (getFragmentsForCompound(c).isEmpty())
					n++;
			set.setResultValue(idx, "compounds w/o hash code", n);
			List<Integer> cmps = new ArrayList<>();
			List<String> cmpsStr = new ArrayList<>();
			for (CFPFragment fragment : fragmentToCompound.keySet())
			{
				cmps.add(fragmentToCompound.get(fragment).size());
				cmpsStr.add(fragmentToCompound.get(fragment).size() + "");
			}
			set.setResultValue(idx, "mean compounds per fragment", DoubleArraySummary.create(cmps));
			CountedSet<String> setStr = CountedSet.create(cmpsStr);
			for (Integer f : new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 })
				set.setResultValue(idx, "fragment with freq " + f + ": ", setStr.getCount(f + ""));

			List<Integer> numFragments = new ArrayList<>();
			for (int c = 0; c < numCompounds; c++)
				numFragments.add(getFragmentsForCompound(c).size());
			set.setResultValue(idx, "mean fragments per compound",
					DoubleArraySummary.create(numFragments));

			if (featureSelection == FeatureSelection.fold)
			{
				estimateCollisions(set, idx, "");
			}
		}
		return set;
	}

	public void estimateCollisions(ResultSet set, int idx, String prefix)
	{
		List<Integer> counts = new ArrayList<Integer>();
		List<String> countsStr = new ArrayList<String>();
		int numCollisions = 0;
		int numBits = 0;
		for (int i = 0; i < hashfoldsize; i++)
			if (collisionMap.containsKey(i))
			{
				if (collisionMap.get(i).size() > 1)
					numCollisions++;
				numBits++;
				counts.add(collisionMap.get(i).size());
				countsStr.add(collisionMap.get(i).size() + "");
			}
			else
				countsStr.add("0");

		//		set.setResultValue(idx, prefix + "collisions", numCollisions + "/" + numBits);
		set.setResultValue(idx, prefix + "collisions", (double) numCollisions / numBits);

		//		set.setResultValue(idx, prefix + "collision ratio", numCollisions / (double) numBits);
		DoubleArraySummary occu = DoubleArraySummary.create(counts);
		set.setResultValue(idx, prefix + "bit-load", occu.getMean());

		//		CountedSet<String> occuStr = CountedSet.create(countsStr);
		//		for (Integer f : new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 })
		//			set.setResultValue(idx, prefix + "bits with #" + f + " fragments", occuStr.getCount(f + ""));
	}

	public String toString()
	{
		return getSummary(false).translate().toNiceString();
	}

	//	LinkedHashMap<Pair, Set<Integer>> pairAdjacent = new LinkedHashMap<>();

	public void toCSVFile(String path, List<String> smiles, List<?> endpoints)
	{
		if (smiles.size() != numCompounds)
			throw new IllegalArgumentException();

		StringBuffer s = new StringBuffer();
		s.append("SMILES,endpoint");
		for (CFPFragment f : fragmentToCompound.keySet())
			s.append("," + f);
		s.append("\n");
		for (int c = 0; c < numCompounds; c++)
		{
			s.append(smiles.get(c) + "," + endpoints.get(c));
			for (CFPFragment f : fragmentToCompound.keySet())
				s.append("," + (getFragmentsForCompound(c).contains(f) ? "1" : "0"));
			s.append("\n");
		}
		FileUtil.writeStringToFile(path, s.toString());
	}

	public String getName()
	{
		String suffix = "";
		if (featureSelection == FeatureSelection.fold)
			suffix = "_" + hashfoldsize;
		else if (featureSelection == FeatureSelection.filt)
			suffix = "_" + hashfoldsize;
		return type + "_" + featureSelection + suffix;
	}

}
