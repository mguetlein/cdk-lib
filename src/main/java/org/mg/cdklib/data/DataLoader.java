package org.mg.cdklib.data;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.mg.cdklib.CDKConverter;
import org.mg.javalib.datamining.ResultSet;
import org.mg.javalib.io.KeyValueFileStore;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.CountedSet;
import org.mg.javalib.util.FileUtil;
import org.mg.javalib.util.FileUtil.CSVFile;
import org.mg.javalib.util.ListUtil;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class DataLoader
{
	public String dataFolder = "data";
	private HashMap<String, CDKDataset> datasets = new HashMap<>();
	private boolean resampleDecoys = true;
	public static boolean VERBOSE = false;

	public static class Source
	{
		String citationKey, citationUrl;

		String getCitationName()
		{
			String s[] = citationKey.split("_");
			if (s.length > 2)
				return s[0].toUpperCase() + " (" + ArrayUtil.last(s) + ")";
			else
				return s[0];
		}

		public Source(String citationKey, String citationUrl)
		{
			this.citationKey = citationKey;
			this.citationUrl = citationUrl;
		}
	}

	public static class WebSource extends Source
	{
		String citationTitle;

		public WebSource(String citationKey, String citationUrl, String citationTitle)
		{
			super(citationKey, citationUrl);
			this.citationTitle = citationTitle;
		}

		public String toLatexURLCitation()
		{
			return "@misc{" + citationKey + ",\n" //
					+ "title = {" + citationTitle + "},\n" //
					+ "howpublished = {\\url{" + citationUrl + "}},\n" //
					+ "note = {Accessed: 2015-XX-XX}\n" //
					+ "}";
		}
	}

	private static void addDatasetWeblink(String dataset, String citationKey, String citationUrl,
			String title)
	{
		Source s = new WebSource(citationKey, citationUrl, title);
		if (!citation.containsKey(dataset))
			citation.put(dataset, new LinkedHashSet<Source>());
		citation.get(dataset).add(s);
	}

	private static void addDatasetCitation(String dataset, String citationKey, String citationUrl)
	{
		Source s = new Source(citationKey, citationUrl);
		if (!citation.containsKey(dataset))
			citation.put(dataset, new LinkedHashSet<Source>());
		citation.get(dataset).add(s);
	}

	public static String OTHER_DATASETS = "Other";
	public static String BALANCED_DATASETS = "Balanced";
	public static String VS_DATASETS = "Virtual-Screening";
	private static HashMap<String, String> datasetCategory = new HashMap<>();
	private static HashMap<String, String> datasetSubCategory = new HashMap<>();
	private static HashMap<String, String> sdfDatasets = new HashMap<>();
	private static HashMap<String, String> sdfEndpoints = new HashMap<>();
	private static HashMap<String, HashSet<Source>> citation = new HashMap<>();
	private static HashMap<String, String> datasetActivityDesc = new HashMap<>();

	static
	{
		String n = null;

		for (String e : new String[] { "MultiCellCall", "SingleCellCall", "Rat", "Mouse", "Hamster",
				"Mutagenicity" }) //"Dog_Primates"
		{
			n = "CPDBAS_" + e;
			sdfDatasets.put(n, "CPDBAS_v5d_1547_20Nov2008.sdf");
			sdfEndpoints.put(n, "ActivityOutcome_CPDBAS_" + e);
			addDatasetWeblink(n, "EPA", "http://www.epa.gov/ncct/dsstox/sdf_cpdbas.html",
					"The Carcinogenic Potency Database");
			addDatasetCitation(n, "gold_supplement_1999",
					"http://toxsci.oxfordjournals.org/content/85/2/747.short");
			datasetActivityDesc.put(n, "carcinogenicity");
			datasetCategory.put(n, BALANCED_DATASETS);
			datasetSubCategory.put(n, "CPDBAS");
		}

		n = "NCTRER";
		addDatasetWeblink(n, "EPA", "http://www.epa.gov/ncct/dsstox/sdf_nctrer.html",
				"Estrogen Receptor Binding Database File");
		addDatasetCitation(n, "fang_structure-activity_2001",
				"http://pubs.acs.org/doi/abs/10.1021/tx000208y");
		datasetActivityDesc.put(n, "Estrogen receptor");
		sdfDatasets.put(n, "NCTRER_v4b_232_15Feb2008.sdf");
		sdfEndpoints.put(n, "ActivityOutcome_NCTRER");
		datasetCategory.put(n, BALANCED_DATASETS);
		datasetSubCategory.put(n, n);

		n = "AMES";
		sdfDatasets.put(n, "cas_4337.ob.noH.sdf");
		sdfEndpoints.put(n, "Ames test categorisation");
		addDatasetWeblink(n, "SD-File", "http://www.cheminformatics.org/datasets/bursi",
				"4337 Structures with AMES Categorisation");
		addDatasetCitation(n, "kazius_derivation_2005",
				"http://pubs.acs.org/doi/abs/10.1021/jm040835a");
		datasetActivityDesc.put(n, "ames test mutagenicity");
		datasetCategory.put(n, BALANCED_DATASETS);
		datasetSubCategory.put(n, n);

		for (String e : new String[] { "vegfr2", "hivrt", "cdk2" })
		{
			n = "DUD_" + e;

			addDatasetCitation(n, "huang_benchmarking_2006",
					"http://pubs.acs.org/doi/abs/10.1021/jm0608356");
			addDatasetWeblink(n, "DUD-Directory", "http://dud.docking.org",
					"Directory of Useful Decoys");
			addDatasetCitation(n, "riniker_heterogeneous_2013",
					"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
			addDatasetWeblink(n, "Benchmarking-Platform",
					"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
			if (e.equals("cdk2"))
				datasetActivityDesc.put(n, "cyclin-dependent kinase");
			else if (e.equals("hivrt"))
				datasetActivityDesc.put(n, "HIV reverse transcriptase");
			else if (e.equals("vegfr2"))
				datasetActivityDesc.put(n, "vascular endothelial growth factor receptor");
			else
				throw new IllegalArgumentException();
			datasetCategory.put(n, VS_DATASETS);
			datasetSubCategory.put(n, "DUD");
		}

		for (Integer i : new Integer[] { 51, 11140, 259, 10434, 10378, 12911, 100, 87, 93, 19905,
				61, 12261, 25, 11365, 100579, 11359, 165, 8, 121, 107, 219, 108, 12670, 104, 105,
				114, 13001, 11575, 11534, 52, 10193, 130, 12252, 17045, 10188, 10498, 10280, 11489,
				12209, 10260, 12952, 10980, 90, 36, 65, 43, 11631, 72, 15, 126 })
		{
			n = "ChEMBL_" + i;
			addDatasetCitation(n, "heikamp_large-scale_2011",
					"http://pubs.acs.org/doi/abs/10.1021/ci200199u");
			addDatasetCitation(n, "riniker_heterogeneous_2013",
					"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
			addDatasetWeblink(n, "Benchmarking-Platform",
					"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
			datasetCategory.put(n, VS_DATASETS);
			datasetSubCategory.put(n, "ChEMBL");
		}
		datasetActivityDesc.put("ChEMBL_" + 101174,
				"pituitary adenylate cyclase-activating polypeptide type I receptor");
		datasetActivityDesc.put("ChEMBL_" + 101395, "IgG receptor FcRn large subunit p51");
		datasetActivityDesc.put("ChEMBL_" + 10102, "5-lipoxygenase activating protein");
		datasetActivityDesc.put("ChEMBL_" + 10144, "bone morphogenetic protein 1");
		datasetActivityDesc.put("ChEMBL_" + 12909, "ileal bile acid transporter");
		datasetActivityDesc.put("ChEMBL_" + 20130, "inhibitor of apoptosis protein");
		datasetActivityDesc.put("ChEMBL_" + 275, "retinoid X receptor alpha");
		datasetActivityDesc.put("ChEMBL_" + 11061, "motilin receptor");
		datasetActivityDesc.put("ChEMBL_" + 10056, "DNA-dependent protein kinase");
		datasetActivityDesc.put("ChEMBL_" + 11096, "sodium/hydrogen exchanger 1");
		datasetActivityDesc.put("ChEMBL_" + 10845, "phospholipase D1");
		datasetActivityDesc.put("ChEMBL_" + 11758, "glucagon-like peptide receptor");
		datasetActivityDesc.put("ChEMBL_" + 11402, "furin");
		datasetActivityDesc.put("ChEMBL_" + 12725, "matriptase");
		datasetActivityDesc.put("ChEMBL_" + 101219, "secreted frizzled-related protein 1");
		datasetActivityDesc.put("ChEMBL_" + 176, "Purinergic receptor P2Y12");
		datasetActivityDesc.put("ChEMBL_" + 10087, "deoxycytidine kinase");
		datasetActivityDesc.put("ChEMBL_" + 100098, "serine/threonine-protein kinase WEE1");
		datasetActivityDesc.put("ChEMBL_" + 10624, "serotonin 5a (5-HT5a) receptor");
		datasetActivityDesc.put("ChEMBL_" + 12659, "prostanoid DP receptor");
		datasetActivityDesc.put("ChEMBL_" + 10582, "cytosolic phospholipase A2");
		datasetActivityDesc.put("ChEMBL_" + 100862, "metastin receptor");
		datasetActivityDesc.put("ChEMBL_" + 117, "somatostatin receptor 2");
		datasetActivityDesc.put("ChEMBL_" + 4,
				"voltage-gated T-type calcium channel alpha-1H subunit");
		datasetActivityDesc.put("ChEMBL_" + 11635, "protein kinase C alpha");
		datasetActivityDesc.put("ChEMBL_" + 11242, "Focal adhesion kinase 1");
		datasetActivityDesc.put("ChEMBL_" + 34, "fibronectin receptor beta");
		datasetActivityDesc.put("ChEMBL_" + 100077, "cell division cycle 7-related protein kinase");
		datasetActivityDesc.put("ChEMBL_" + 193, "coagulation factor IX");
		datasetActivityDesc.put("ChEMBL_" + 80, "FK506-binding protein 1A");
		datasetActivityDesc.put("ChEMBL_" + 165, "HERG");
		datasetActivityDesc.put("ChEMBL_" + 10193, "carbonic anhydrase I");
		datasetActivityDesc.put("ChEMBL_" + 15, "carbonic anhydrase II");
		datasetActivityDesc.put("ChEMBL_" + 11489, "11-beta-hydroxysteroid dehydrogenase 1");
		datasetActivityDesc.put("ChEMBL_" + 121, "serotonin transporter");
		datasetActivityDesc.put("ChEMBL_" + 72, "dopamine D2 receptor");
		datasetActivityDesc.put("ChEMBL_" + 259, "cannabinoid CB2 receptor");
		datasetActivityDesc.put("ChEMBL_" + 10188, "MAP kinase p38 alpha");
		datasetActivityDesc.put("ChEMBL_" + 108, "serotonin 2c (5-HT2c) receptor");
		datasetActivityDesc.put("ChEMBL_" + 12952, "carbonic anhydrase IX");
		datasetActivityDesc.put("ChEMBL_" + 93, "acetylcholinesterase");
		datasetActivityDesc.put("ChEMBL_" + 10980, "vascular endothelial growth factor receptor 2");
		datasetActivityDesc.put("ChEMBL_" + 19905, "melanin-concentrating hormone receptor 1");
		datasetActivityDesc.put("ChEMBL_" + 107, "serotonin 2a (5-HT2a) receptor");
		datasetActivityDesc.put("ChEMBL_" + 87, "cannabinoid CB1 receptor");
		datasetActivityDesc.put("ChEMBL_" + 17045, "cytochrome P450 3A4");
		datasetActivityDesc.put("ChEMBL_" + 11140, "dipeptidyl peptidase IV");
		datasetActivityDesc.put("ChEMBL_" + 114, "adenosine A1 receptor");
		datasetActivityDesc.put("ChEMBL_" + 90, "dopamine D4 receptor");
		datasetActivityDesc.put("ChEMBL_" + 100166, "kinesin-like protein 1");
		datasetActivityDesc.put("ChEMBL_" + 13001, "matrix metalloproteinase-2");
		datasetActivityDesc.put("ChEMBL_" + 104, "monoamine oxidase B");
		datasetActivityDesc.put("ChEMBL_" + 65, "cytochrome P450 19A1");
		datasetActivityDesc.put("ChEMBL_" + 61, "muscarinic acetylcholine receptor M1");
		datasetActivityDesc.put("ChEMBL_" + 10280, "histamine H3 receptor");
		datasetActivityDesc.put("ChEMBL_" + 51, "serotonin 1a (5-HT1a) receptor");
		datasetActivityDesc.put("ChEMBL_" + 100, "norepinephrine transporter");
		datasetActivityDesc.put("ChEMBL_" + 10260, "vanilloid receptor");
		datasetActivityDesc.put("ChEMBL_" + 52, "alpha-2a adrenergic receptor");
		datasetActivityDesc.put("ChEMBL_" + 11365, "cytochrome P450 2D6");
		datasetActivityDesc.put("ChEMBL_" + 11359, "phosphodiesterase 4D");
		datasetActivityDesc.put("ChEMBL_" + 28, "thymidylate synthase");
		datasetActivityDesc.put("ChEMBL_" + 11536, "ghrelin receptor");
		datasetActivityDesc.put("ChEMBL_" + 8, "tyrosine-protein kinase ABL");
		datasetActivityDesc.put("ChEMBL_" + 10434, "tyrosine-protein kinase SRC");
		datasetActivityDesc.put("ChEMBL_" + 12670, "tyrosine-protein kinase receptor FLT3");
		datasetActivityDesc.put("ChEMBL_" + 20014, "serine/threonine-protein kinase Aurora-A");
		datasetActivityDesc.put("ChEMBL_" + 234, "insulin-like growth factor I receptor");
		datasetActivityDesc.put("ChEMBL_" + 12261, "c-Jun N-terminal kinase 1");
		datasetActivityDesc.put("ChEMBL_" + 12209, "carbonic anhydrase XII");
		datasetActivityDesc.put("ChEMBL_" + 25, "glucocorticoid receptor");
		datasetActivityDesc.put("ChEMBL_" + 36, "progesterone receptor");
		datasetActivityDesc.put("ChEMBL_" + 43, "beta-2 adrenergic receptor");
		datasetActivityDesc.put("ChEMBL_" + 219, "muscarinic acetylcholine receptor M3");
		datasetActivityDesc.put("ChEMBL_" + 130, "dopamine D3 receptor");
		datasetActivityDesc.put("ChEMBL_" + 105, "serotonin 1d (5-HT1d) receptor");
		datasetActivityDesc.put("ChEMBL_" + 11336, "neuropeptide Y receptor type 5");
		datasetActivityDesc.put("ChEMBL_" + 20174, "G protein-coupled receptor");
		datasetActivityDesc.put("ChEMBL_" + 126, "cyclooxygenase-2");
		datasetActivityDesc.put("ChEMBL_" + 11225, "renin");
		datasetActivityDesc.put("ChEMBL_" + 12252, "beta-secretase 1");
		datasetActivityDesc.put("ChEMBL_" + 11682, "glycine transporter 1");
		datasetActivityDesc.put("ChEMBL_" + 134, "vasopressin V1a receptor");
		datasetActivityDesc.put("ChEMBL_" + 116, "oxytocin receptor");
		datasetActivityDesc.put("ChEMBL_" + 11265, "somatostatin receptor 5");
		datasetActivityDesc.put("ChEMBL_" + 10475, "neuropeptide Y receptor type 1");
		datasetActivityDesc.put("ChEMBL_" + 12679, "C5a anaphylatoxin chemotactic receptor");
		datasetActivityDesc.put("ChEMBL_" + 10579, "C-C chemokine receptor type 4");
		datasetActivityDesc.put("ChEMBL_" + 11575, "C-C chemokine receptor type 2");
		datasetActivityDesc.put("ChEMBL_" + 18061, "sodium channel protein type IX alpha subunit");
		datasetActivityDesc.put("ChEMBL_" + 237, "leukotriene A4 hydrolase");
		datasetActivityDesc.put("ChEMBL_" + 276, "phosphodiesterase 4A");
		datasetActivityDesc.put("ChEMBL_" + 11534, "cathepsin S");
		datasetActivityDesc.put("ChEMBL_" + 10198, "voltage-gated potassium channel subunit Kv1.5");
		datasetActivityDesc.put("ChEMBL_" + 10498, "cathepsin L");
		datasetActivityDesc.put("ChEMBL_" + 12911, "cytochrome P450 2C9");
		datasetActivityDesc.put("ChEMBL_" + 12968, "orexin receptor 2");
		datasetActivityDesc.put("ChEMBL_" + 100579, "nicotinic acid receptor 1");
		datasetActivityDesc.put("ChEMBL_" + 100126, "serine/threonine-protein kinase B-raf");
		datasetActivityDesc.put("ChEMBL_" + 10378, "cathepsin B");
		datasetActivityDesc.put("ChEMBL_" + 10417, "P2X purinoceptor 7");
		datasetActivityDesc.put("ChEMBL_" + 10752,
				"inhibitor of nuclear factor kappa B kinase beta subunit");
		datasetActivityDesc.put("ChEMBL_" + 10773, "interleukin-8 receptor B");
		datasetActivityDesc.put("ChEMBL_" + 11631, "sphingosine 1-phosphate receptor Edg-1");
		datasetActivityDesc.put("ChEMBL_" + 10927, "urotensin II receptor");
		datasetActivityDesc.put("ChEMBL_" + 11085, "melatonin receptor 1B");
		datasetActivityDesc.put("ChEMBL_" + 11442, "liver glycogen phosphorylase");
		datasetActivityDesc.put("ChEMBL_" + 11279, "metabotropic glutamate receptor 1");
		datasetActivityDesc.put("ChEMBL_" + 11488, "estradiol 17-beta-dehydrogenase 3");
		datasetActivityDesc.put("ChEMBL_" + 12840, "macrophage colony stimulating factor receptor");

		for (Integer i : new Integer[] { 644, 713, 859, 737, 852, 692, 652, 689, 712, 600, 466, 832,
				858, 810, 733, 548 })
		{
			n = "MUV_" + i;
			addDatasetCitation(n, "rohrer_maximum_2009",
					"http://pubs.acs.org/doi/abs/10.1021/ci8002649");
			addDatasetWeblink(n, "MUV-Datasets",
					"http://www.pharmchem.tu-bs.de/lehre/baumann/MUV.html",
					"Maximum Unbiased Validation (MUV) Datasets for Virtual Screening");
			addDatasetCitation(n, "riniker_heterogeneous_2013",
					"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
			addDatasetWeblink(n, "Benchmarking-Platform",
					"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
			datasetCategory.put(n, VS_DATASETS);
			datasetSubCategory.put(n, "MUV");
		}
		datasetActivityDesc.put("MUV_" + 466, "S1P1 rec. (GPCR) Agonist");
		datasetActivityDesc.put("MUV_" + 548, "PKA (Kinase) Inhibitor");
		datasetActivityDesc.put("MUV_" + 600, "SF1 (Nuclear Receptor) Inhibitor");
		datasetActivityDesc.put("MUV_" + 644, "Rho-Kinase2 Inhibitor");
		datasetActivityDesc.put("MUV_" + 652, "HIV RT-RNase Inhibitor");
		datasetActivityDesc.put("MUV_" + 689, "Eph rec. A4 (Rec. Tyr. Kinase) Inhibitor");
		datasetActivityDesc.put("MUV_" + 692, "SF1 (Nuclear Receptor) Agonist");
		datasetActivityDesc.put("MUV_" + 712, "HSP 90 (Chaperone) Inhibitor");
		datasetActivityDesc.put("MUV_" + 713, "ER-a-Coact. Bind. (PPI) Inhibitor");
		datasetActivityDesc.put("MUV_" + 733, "ER-ß-Coact. Bind. (PPI) Inhibitor");
		datasetActivityDesc.put("MUV_" + 737, "ER-a-Coact. Bind. (PPI) Potentiator");
		datasetActivityDesc.put("MUV_" + 810, "FAK (Kinase) Inhibitor");
		datasetActivityDesc.put("MUV_" + 832, "Cathepsin G (Protease) Inhibitor");
		datasetActivityDesc.put("MUV_" + 846, "FXIa (Protease) Inhibitor");
		datasetActivityDesc.put("MUV_" + 852, "FXIIa (Protease) Inhibitor");
		datasetActivityDesc.put("MUV_" + 858, "D1 rec. (GPCR) Allosteric Modulator");
		datasetActivityDesc.put("MUV_" + 859, "M1 rec. (GPCR) Allosteric Modulator");

		//		n = "REID";
		//		datasetActivityDesc.put(n, "norm induction >= 11");
		//		sdfDatasets.put(n, "all_isoxazoles.sdf");
		//		sdfEndpoints.put(n, "activity");
		//		addDatasetCitation(n, "-", "-");
		//		addDatasetWeblink(n, "-", "-", "-");
		//		datasetCategory.put(n, OTHER_DATASETS);
		//		datasetSubCategory.put(n, n);

		for (int i : new int[] { 2, 3, 4, 11 })
		{
			n = "REID-" + i;
			datasetActivityDesc.put(n, "norm induction >= " + i);
			sdfDatasets.put(n, "reid-6.sdf");
			sdfEndpoints.put(n, "activity" + i);
			addDatasetCitation(n, "-", "-");
			addDatasetWeblink(n, "-", "-", "-");
			datasetCategory.put(n, OTHER_DATASETS);
			datasetSubCategory.put(n, n);
		}
		n = "LTKB";
		datasetActivityDesc.put(n, "drug-induced liver injury (DILI)");
		addDatasetCitation(n, "-", "-");
		addDatasetWeblink(n, "Liver Toxicity Knowledge Base (LTKB)",
				"http://www.fda.gov/ScienceResearch/BioinformaticsTools/LiverToxicityKnowledgeBase/ucm2024036.htm",
				"Liver Toxicity Knowledge Base (LTKB)");
		datasetCategory.put(n, OTHER_DATASETS);
		datasetSubCategory.put(n, n);
	}

	public String[] reducedDatasets()
	{
		return FileUtil.readStringFromFile(dataFolder + "/endpoints_reduced.txt").split("\n");
	}

	public String[] allDatasets()
	{
		return FileUtil.readStringFromFile(dataFolder + "/endpoints.txt").split("\n");
	}

	public String[] allDatasetsSorted()
	{
		String[] s = FileUtil.readStringFromFile(dataFolder + "/endpoints.txt").split("\n");
		Arrays.sort(s, CFPDataComparator);
		return s;
	}

	public String[] cpdbDatasets()
	{
		return filteredDatasets("CPDB");
	}

	public String[] chemblDatasets()
	{
		return filteredDatasets("ChEMBL");
	}

	public String[] muvDatasets()
	{
		return filteredDatasets("MUV");
	}

	public String[] dudDatasets()
	{
		return filteredDatasets("DUD");
	}

	public String[] filteredDatasets(final String prefix)
	{
		String[] s = FileUtil.readStringFromFile(dataFolder + "/endpoints.txt").split("\n");
		return ArrayUtil.filter(String.class, s, new ListUtil.Filter<String>()
		{
			@Override
			public boolean accept(String s)
			{
				return s.startsWith(prefix);
			}
		});
	}

	public String[] balancedDatasets()
	{
		String[] s = FileUtil.readStringFromFile(dataFolder + "/endpoints.txt").split("\n");
		return ArrayUtil.filter(String.class, s, new ListUtil.Filter<String>()
		{
			@Override
			public boolean accept(String s)
			{
				return datasetCategory.get(s) == BALANCED_DATASETS;
			}
		});
	}

	public static final DataLoader INSTANCE = new DataLoader(
			System.getProperty("user.home") + "/results/cdklib/data");

	private DataLoader(String dataFolder)
	{
		this.dataFolder = dataFolder;
	}

	//	public boolean exists(String name)
	//	{
	//		if (datasets.containsKey(name))
	//			return true;
	//		if (sdfDatasets.containsKey(name))
	//			return new File(dataFolder + File.separator + sdfDatasets.get(name)).exists();
	//		return new File(dataFolder + File.separator + name + ".csv").exists();
	//	}

	public ResultSet getInfo(boolean cite, String... ids)
	{
		ResultSet set = new ResultSet();
		set.setNicePropery("size", "compounds");
		for (String n : ids)
		{
			//			System.out.println(n);

			int rIdx = set.addResult();
			set.setResultValue(rIdx, "category", datasetCategory.get(n));
			set.setResultValue(rIdx, "name", n);//.replaceAll("_", " "));
			CDKDataset d = getDataset(n);
			String classV[] = getClassValues(d.endpoints);
			int activeIdx = getActiveIdx(classV);
			set.setResultValue(rIdx, "size", d.getSmiles().size());
			CountedSet<String> endp = CountedSet.create(d.getEndpoints());
			set.setResultValue(rIdx, "active", endp.getCount(classV[activeIdx]));
			set.setResultValue(rIdx, "in-active", endp.getCount(classV[1 - activeIdx]));
			//			set.setResultValue(rIdx, "#in-active", endp.getCount(classV[1 - activeIdx]));
			//			set.setResultValue(rIdx, "activity", ArrayUtil.toString(classV));
			set.setResultValue(rIdx, "target", datasetActivityDesc.get(n));

			//			set.setResultValue(rIdx, "dataset-weblink", CollectionUtil.toString(datasetWeblinks.get(n)));

			String source;
			if (cite)
			{
				String cit = "\\cite{";
				for (Source s : citation.get(n))
					if (!(s instanceof WebSource))
						cit += s.citationKey + ",";
				cit = cit.substring(0, cit.length() - 1) + "}";
				source = cit;
			}
			else
			{
				source = "";
				for (Source s : citation.get(n))
					if (s instanceof WebSource)
						source += "\\url{" + s.citationUrl + "}, ";
				source = source.substring(0, source.length() - 2);
			}
			set.setResultValue(rIdx, "source", source);
		}
		return set;
	}

	public ResultSet getCategoryInfo(boolean cite, String... name)
	{
		ResultSet set = getInfo(cite, name);

		for (int idx = 0; idx < set.getNumResults(); idx++)
		{
			String n = set.getResultValue(idx, "name").toString();
			set.setResultValue(idx, "category", datasetCategory.get(n));
			set.setResultValue(idx, "subCategory", datasetSubCategory.get(n));
			set.setResultValue(idx, "numDatasets", "1");
		}

		List<String> props = new ArrayList<>(set.getProperties());
		props.add(1, props.remove(props.size() - 1));
		props.add(1, props.remove(props.size() - 1));
		set.sortProperties(props);
		set = set.join(new String[] { "category", "subCategory", "source" },
				new String[] { "name", "target" }, null);

		//		System.out.println(set.toNiceString());

		for (int i = 0; i < set.getNumResults(); i++)
		{
			set.setResultValue(i, "numDatasets",
					set.getResultValue(i, "numDatasets").toString().split("/").length);
			for (String p : new String[] { "size", "active", "in-active" })
				if (((Integer) set.getResultValue(i, "numDatasets")) > 1)
				{
					String mean = new DecimalFormat("#.#")
							.format((Double) set.getResultValue(i, p));
					//					if (!mean.contains("."))
					//						mean += ".0";
					set.setResultValue(i, p, mean);
				}
			if (set.getResultValue(i, "category").equals(set.getResultValue(i, "subCategory")))
				set.setResultValue(i, "subCategory", "-");
			//				set.setResultValue(i, "category",
			//						set.getResultValue(i, "category") + " -- " + set.getResultValue(i, "subCategory"));
		}
		//		set.removePropery("subCategory");
		set.setNicePropery("category", "type");
		set.setNicePropery("subCategory", "dataset/group");
		set.setNicePropery("numDatasets", "num");
		set.setNicePropery("size", "compounds");
		//		System.out.println(set.toNiceString());

		if (!cite)
		{
			set.removePropery("category");
			set.removePropery("numDatasets");
			set.removePropery("size");
			set.removePropery("active");
			set.removePropery("in-active");
		}

		return set;
	}

	public void setResampleDecoys(boolean resampleDecoys)
	{
		this.resampleDecoys = resampleDecoys;
	}

	@Deprecated
	public CDKDataset getDataset(String name, int run)
	{
		return getDataset(name);
	}

	KeyValueFileStore<String, CDKDataset> parsedDatasets;

	public CDKDataset getDataset(String name)
	{
		if (!datasets.containsKey(name))
		{
			try
			{
				if (parsedDatasets == null)
					parsedDatasets = new KeyValueFileStore<>(dataFolder + "/parsed", false, true,
							null, true);
				CDKDataset data;
				if (!parsedDatasets.contains(name))
				{
					if (sdfDatasets.containsKey(name))
						data = getDatasetFromSDF(name);
					else if (new File(dataFolder + File.separator + name + ".csv").exists())
						data = getDatasetFromCSV(name);
					else if (new File(dataFolder + File.separator + name + ".smi").exists())
						data = getDatasetFromSMILES(name);
					else
						throw new RuntimeException(
								"dataset not found: " + name + " in " + dataFolder);
					parsedDatasets.store(name, data);
				}
				else
					data = parsedDatasets.get(name);
				datasets.put(name, data);
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}
		}
		return datasets.get(name);
	}

	public boolean hasDuplicates(String name)
	{
		Integer run = 1;
		if (sdfDatasets.containsKey(name))
			run = null;
		if (run != null)
		{
			if (!resampleDecoys)
				run = 1;
			name = name + "_r" + String.format("%02d", run);
		}

		try
		{
			if (sdfDatasets.containsKey(name))
				return hasDuplicatesSDF(name);
			else
				return hasDuplicatesCSV(name);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}

	public static CDKDataset getDatasetFromSDF(String name, String path)
			throws CDKException, FileNotFoundException, IOException
	{
		return getDatasetFromSDF(name, path, null);
	}

	public static CDKDataset getDatasetFromSDF(String name, String path, String endpoint)
			throws CDKException, FileNotFoundException, IOException
	{
		if (VERBOSE)
			System.err.println("data-loader> " + endpoint);

		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();

		ISimpleChemObjectReader reader = new ReaderFactory()
				.createReader(new InputStreamReader(new FileInputStream(path)));
		IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
		int invalidCompound = 0;
		int missingEndpoint = 0;
		for (IAtomContainer a : ChemFileManipulator.getAllAtomContainers(content))
		{
			if (a.getAtomCount() == 0)
				invalidCompound++;
			else if (endpoint != null && (a.getProperty(endpoint) == null
					|| a.getProperty(endpoint).toString().equals("unspecified")
					|| a.getProperty(endpoint).toString().equals("blank")
					|| a.getProperty(endpoint).toString().equals("inconclusive")))
				missingEndpoint++;
			else
			{
				String smi = new SmilesGenerator().create(a);
				//				CDKUtil.setMolForSmiles(smi, a);

				//				try
				//				{
				//					IAtomContainer m2 = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smi);
				//					if (m2.getAtomCount() == 0)
				//						throw new RuntimeException("num atoms 0");
				//					if (m2.getAtomCount() != a.getAtomCount())
				//						throw new RuntimeException("num atoms " + a.getAtomCount() + " != " + m2.getAtomCount());
				//				}
				//				catch (Exception e)
				//				{
				//					e.printStackTrace();
				//					System.exit(1);
				//				}
				smiles.add(smi);
				if (endpoint != null)
					endpoints.add(a.getProperty(endpoint).toString());
				else
					endpoints.add("n/a");
			}
		}
		reader.close();

		List<String> warnings = new ArrayList<>();
		if (invalidCompound > 0)
			warnings.add("Removed " + invalidCompound
					+ " compounds that could not be read by the CDK library.");
		if (missingEndpoint > 0)
			warnings.add("Removed " + missingEndpoint
					+ " compounds with missing/invalid enpoint values.");

		return createDataset(name, smiles, endpoints, warnings);
	}

	private CDKDataset getDatasetFromSDF(String name)
			throws FileNotFoundException, CDKException, IOException
	{
		return getDatasetFromSDF(name, dataFolder + File.separator + sdfDatasets.get(name),
				sdfEndpoints.get(name));
	}

	private static CDKDataset createDataset(String name, List<String> smiles,
			List<String> endpoints, List<String> warnings) throws CDKException
	{
		HashMap<String, HashSet<String>> uniqToActivity = new HashMap<>();
		int idx = 0;
		for (String smi : smiles)
		{
			//			System.err.println("smiles: " + smi);
			String uniq = CDKConverter.toAbsoluteSmiles(smi);
			if (!uniqToActivity.containsKey(uniq))
				uniqToActivity.put(uniq, new HashSet<String>());
			uniqToActivity.get(uniq).add(endpoints.get(idx));
			idx++;
		}

		List<String> uSmiles = new ArrayList<>();
		List<String> uEndpoints = new ArrayList<>();
		HashMap<String, Integer> uniqAdded = new HashMap<>();
		int skipEq = 0;
		int skipDiff = 0;
		for (int i = 0; i < smiles.size(); i++)
		{
			//					System.out.println(i + " " + smiles.get(i));
			String uniq = CDKConverter.toAbsoluteSmiles(smiles.get(i));
			if (uniqAdded.containsKey(uniq))
			{
				skipEq++;
				if (VERBOSE)
					System.err.println("createDataset> skip equal value " + uniqToActivity.get(uniq)
							+ " duplicate " + uniq + ":\n" + i + " " + smiles.get(i)
							+ "\nis equal to:\n" + uniqAdded.get(uniq) + " "
							+ smiles.get(uniqAdded.get(uniq)));
			}
			else if (uniqToActivity.get(uniq).size() > 1)
			{
				skipDiff++;
				if (VERBOSE)
					System.err.println("createDataset> skip different value duplicate: " + i + " "
							+ uniqToActivity.get(uniq) + " " + smiles.get(i) + " " + uniq);
			}
			else
			{
				uSmiles.add(smiles.get(i));
				uEndpoints.add(endpoints.get(i));
				uniqAdded.put(uniq, i);
			}
		}
		if (skipEq > 0)
			warnings.add("Removed " + skipEq
					+ " duplicate occurences of compounds (with equal endpoint values).");
		if (skipDiff > 0)
			warnings.add("Removed " + skipDiff
					+ " compounds that occured multiple times with different endpoint values.");

		return new CDKDataset(name, uSmiles, uEndpoints, warnings);
	}

	private boolean hasDuplicatesSDF(String name) throws Exception
	{
		List<String> smiles = new ArrayList<>();
		ISimpleChemObjectReader reader = new ReaderFactory().createReader(new InputStreamReader(
				new FileInputStream(dataFolder + File.separator + sdfDatasets.get(name))));
		IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
		String endpoint = sdfEndpoints.get(name);
		for (IAtomContainer a : ChemFileManipulator.getAllAtomContainers(content))
			if (a.getAtomCount() > 0 && a.getProperty(endpoint) != null
					&& !a.getProperty(endpoint).toString().equals("unspecified")
					&& !a.getProperty(endpoint).toString().equals("blank")
					&& !a.getProperty(endpoint).toString().equals("inconclusive"))
			{
				String smi = new SmilesGenerator().create(a);
				smiles.add(CDKConverter.toAbsoluteSmiles(smi));
			}
		reader.close();
		return CountedSet.create(smiles).getMaxCount(false) > 1;
	}

	private CDKDataset getDatasetFromCSV(String name) throws CDKException
	{
		return getDatasetFromCSV(name, dataFolder + File.separator + name + ".csv", 0, 1);
	}

	public static CDKDataset getDatasetFromCSV(String name, String path, int smilesCol,
			int endpointsCol) throws CDKException
	{
		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();

		CSVFile csv = FileUtil.readCSV(path);
		for (int i = 1; i < csv.content.size(); i++)
		{
			smiles.add(csv.content.get(i)[smilesCol]);
			endpoints.add(csv.content.get(i)[endpointsCol]);
		}

		return createDataset(name, smiles, endpoints, new ArrayList<String>());
	}

	private CDKDataset getDatasetFromSMILES(String name) throws CDKException
	{
		return getDatasetFromSMILES(name, dataFolder + File.separator + name + ".smi");
	}

	public static CDKDataset getDatasetFromSMILES(String name, String path) throws CDKException
	{
		return getDatasetFromSMILES(name, path, false);
	}

	public static CDKDataset getDatasetFromSMILES(String name, String path, boolean noEndpoints)
			throws CDKException
	{
		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();
		int numInvalid = 0;

		for (String s : FileUtil.readStringFromFile(path).split("\n"))
		{
			String ss[] = s.split("\\s");
			try
			{
				CDKConverter.toAbsoluteSmiles(ss[0]);
				smiles.add(ss[0]);
				endpoints.add((ss.length > 1 && !noEndpoints) ? ss[1] : "");
			}
			catch (CDKException e)
			{
				numInvalid++;
				if (VERBOSE)
					System.err.println("createDataset> could not read smiles: " + ss[0]);
			}
		}

		List<String> warnings = new ArrayList<>();
		if (numInvalid > 0)
			warnings.add("Could not read " + numInvalid + " smiles.");
		return createDataset(name, smiles, endpoints, new ArrayList<String>());
	}

	private boolean hasDuplicatesCSV(String name) throws CDKException
	{
		List<String> smiles = new ArrayList<>();
		CSVFile csv = FileUtil.readCSV(dataFolder + File.separator + name + ".csv");
		for (int i = 1; i < csv.content.size(); i++)
			smiles.add(CDKConverter.toAbsoluteSmiles(csv.content.get(i)[0]));
		return CountedSet.create(smiles).getMaxCount(false) > 1;
	}

	public static Comparator<Object> CFPDataComparator = new Comparator<Object>()
	{
		@Override
		public int compare(Object o1, Object o2)
		{
			String s1 = o1.toString();
			String s2 = o2.toString();
			if (s1.equals("NCTRER"))
				if (s2.startsWith("ChEMBL") || s2.startsWith("MUV") || s2.startsWith("DUD"))
					return -1;
				else
					return 1;
			else if (s2.equals("NCTRER"))
				if (s1.startsWith("ChEMBL") || s1.startsWith("MUV") || s1.startsWith("DUD"))
					return 1;
				else
					return -1;
			else if (s1.startsWith("ChEMBL") && s2.startsWith("ChEMBL"))
			{
				Integer i1 = Integer.parseInt(s1.substring(s1.indexOf("_") + 1));
				Integer i2 = Integer.parseInt(s2.substring(s2.indexOf("_") + 1));
				return i1.compareTo(i2);
			}
			else
				return s1.compareTo(s2);
		}
	};

	public static void main(String[] args)
	{
		DataLoader d = INSTANCE;
		//		System.out.println(d.getCategoryInfo(d.allDatasets()).toLatexTable());

		//		d.getInfo("NCTRER", "CPDBAS_Rat", "ChEMBL_61", "DUD_vegfr2", "DUD_hivrt", "DUD_cdk2", "MUV_644", "MUV_713",
		//				"MUV_859", "AMES");
		String data[] = d.allDatasets();
		//		Arrays.sort(data, CFPDataComparator);
		//		System.out.println(ArrayUtil.toString(data));
		data = new String[] { "CPDBAS_Mouse" };
		for (String dat : data)
		{
			CDKDataset da = d.getDataset(dat);
			if (!da.warnings.isEmpty())
			{
				System.out.println(dat + " Warnings:");
				for (String warn : da.warnings)
					System.out.println("* " + warn);
				System.out.println();
			}
		}
	}

	public static String[] getClassValues(List<String> endpoints)
	{
		String classValues[] = ArrayUtil.removeDuplicates(ArrayUtil.toArray(endpoints));
		Arrays.sort(classValues);
		return classValues;
	}

	public static Integer getActiveIdx(String[] classValues)
	{
		Integer activeIdx = null;
		for (int i = 0; i < classValues.length; i++)
			if (classValues[i].equals("active") || classValues[i].equals("mutagen")
					|| classValues[i].equals("1") || classValues[i].equals("most-concern"))
				activeIdx = i;
		if (activeIdx == null)
			throw new IllegalStateException("what is active? " + ArrayUtil.toString(classValues));
		return activeIdx;
	}

	public String getDatasetEndpoint(String name)
	{
		return datasetActivityDesc.get(name);
	}

	public Set<String> getDatasetURLs(String name)
	{
		HashSet<String> urls = new LinkedHashSet<>();
		for (Source s : citation.get(name))
			urls.add(s.citationUrl);
		return urls;
	}

	public Map<String, String> getModelDatasetCitations(String modelId)
	{
		HashMap<String, String> map = new LinkedHashMap<>();
		for (Source s : citation.get(modelId))
			map.put(s.getCitationName(), s.citationUrl);
		return map;
	}

	public static Set<String> listCategoryDatasets(String category)
	{
		HashSet<String> datasets = new HashSet<>();
		for (String s : datasetCategory.keySet())
			if (datasetCategory.get(s).equals(category))
				datasets.add(s);
		return datasets;
	}

	public static Set<String> listSubCategoryDatasets(String category)
	{
		HashSet<String> datasets = new HashSet<>();
		for (String s : datasetSubCategory.keySet())
			if (datasetSubCategory.get(s).equals(category))
				datasets.add(s);
		return datasets;
	}

	public boolean isVirtualScreeningDataset(String name)
	{
		return datasetCategory.get(name) == VS_DATASETS;
	}

}
