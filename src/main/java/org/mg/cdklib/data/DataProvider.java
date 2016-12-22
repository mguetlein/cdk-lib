package org.mg.cdklib.data;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.mg.cdklib.CDKConverter;
import org.mg.javalib.datamining.ResultSet;
import org.mg.javalib.io.KeyValueFileStore;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.CountedSet;
import org.mg.javalib.util.FileUtil;
import org.mg.javalib.util.FileUtil.CSVFile;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class DataProvider
{
	public static enum DataID
	{
		AMES, NCTRER, // 
		CPDBAS_Mouse, CPDBAS_MultiCellCall, CPDBAS_Mutagenicity, CPDBAS_Rat, CPDBAS_SingleCellCall,
		CPDBAS_Hamster, CPDBAS_Dog_Primates, //
		DUD_cdk2, DUD_hivrt, DUD_vegfr2, // 
		ChEMBL_100, ChEMBL_100579, ChEMBL_10188, ChEMBL_10193, ChEMBL_10260, ChEMBL_10280,
		ChEMBL_10378, ChEMBL_104, ChEMBL_10434, ChEMBL_10498, ChEMBL_105, ChEMBL_107, ChEMBL_108,
		ChEMBL_10980, ChEMBL_11140, ChEMBL_11359, ChEMBL_11365, ChEMBL_114, ChEMBL_11489,
		ChEMBL_11534, ChEMBL_11575, ChEMBL_11631, ChEMBL_121, ChEMBL_12209, ChEMBL_12252,
		ChEMBL_12261, ChEMBL_126, ChEMBL_12670, ChEMBL_12911, ChEMBL_12952, ChEMBL_130,
		ChEMBL_13001, ChEMBL_15, ChEMBL_165, ChEMBL_17045, ChEMBL_19905, ChEMBL_219, ChEMBL_25,
		ChEMBL_259, ChEMBL_36, ChEMBL_43, ChEMBL_51, ChEMBL_52, ChEMBL_61, ChEMBL_65, ChEMBL_72,
		ChEMBL_8, ChEMBL_87, ChEMBL_90, ChEMBL_93, //
		MUV_466, MUV_548, MUV_600, MUV_644, MUV_652, MUV_689, MUV_692, MUV_712, MUV_713, MUV_733,
		MUV_737, MUV_810, MUV_832, MUV_852, MUV_858, MUV_859, //
		LTKB;
	}

	public static String DATA_FOLDER = System.getProperty("user.home") + "/results/cdklib/data";

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

	public static String OTHER_DATASETS = "Other";
	public static String BALANCED_DATASETS = "Balanced";
	public static String VS_DATASETS = "Virtual-Screening";
	private static HashMap<DataID, CDKDataset> datasets = new HashMap<>();
	private static Set<DataID> cfpDatasets = new HashSet<>();
	private static HashMap<DataID, String> datasetCategory = new HashMap<>();
	private static HashMap<DataID, String> datasetSubCategory = new HashMap<>();
	private static HashMap<DataID, String> sdfDatasets = new HashMap<>();
	private static HashMap<DataID, String> sdfEndpoints = new HashMap<>();
	private static HashMap<DataID, HashSet<Source>> citation = new HashMap<>();
	private static HashMap<DataID, String> datasetActivityDesc = new HashMap<>();

	private static void addDatasetWeblink(DataID dataset, String citationKey, String citationUrl,
			String title)
	{
		Source s = new WebSource(citationKey, citationUrl, title);
		if (!citation.containsKey(dataset))
			citation.put(dataset, new LinkedHashSet<Source>());
		citation.get(dataset).add(s);
	}

	private static void addDatasetCitation(DataID dataset, String citationKey, String citationUrl)
	{
		Source s = new Source(citationKey, citationUrl);
		if (!citation.containsKey(dataset))
			citation.put(dataset, new LinkedHashSet<Source>());
		citation.get(dataset).add(s);
	}

	static
	{
		for (DataID n : DataID.values())
		{
			if (n.toString().startsWith("CPDBAS_"))
			{
				sdfDatasets.put(n, "CPDBAS_v5d_1547_20Nov2008.sdf");
				sdfEndpoints.put(n, "ActivityOutcome_" + n.toString());
				addDatasetWeblink(n, "EPA", "http://www.epa.gov/ncct/dsstox/sdf_cpdbas.html",
						"The Carcinogenic Potency Database");
				addDatasetCitation(n, "gold_supplement_1999",
						"http://toxsci.oxfordjournals.org/content/85/2/747.short");
				datasetActivityDesc.put(n, "carcinogenicity");
				// two small
				if (n != DataID.CPDBAS_Hamster && n != DataID.CPDBAS_Dog_Primates)
					cfpDatasets.add(n);
				datasetCategory.put(n, BALANCED_DATASETS);
				datasetSubCategory.put(n, "CPDBAS");
			}
			else if (n == DataID.NCTRER)
			{
				n = DataID.valueOf("NCTRER");
				addDatasetWeblink(n, "EPA", "http://www.epa.gov/ncct/dsstox/sdf_nctrer.html",
						"Estrogen Receptor Binding Database File");
				addDatasetCitation(n, "fang_structure-activity_2001",
						"http://pubs.acs.org/doi/abs/10.1021/tx000208y");
				datasetActivityDesc.put(n, "Estrogen receptor");
				sdfDatasets.put(n, "NCTRER_v4b_232_15Feb2008.sdf");
				sdfEndpoints.put(n, "ActivityOutcome_NCTRER");
				cfpDatasets.add(n);
				datasetCategory.put(n, BALANCED_DATASETS);
				datasetSubCategory.put(n, n.toString());
			}
			else if (n == DataID.AMES)
			{
				sdfDatasets.put(n, "cas_4337.ob.noH.sdf");
				sdfEndpoints.put(n, "Ames test categorisation");
				addDatasetWeblink(n, "SD-File", "http://www.cheminformatics.org/datasets/bursi",
						"4337 Structures with AMES Categorisation");
				addDatasetCitation(n, "kazius_derivation_2005",
						"http://pubs.acs.org/doi/abs/10.1021/jm040835a");
				datasetActivityDesc.put(n, "ames test mutagenicity");
				cfpDatasets.add(n);
				datasetCategory.put(n, BALANCED_DATASETS);
				datasetSubCategory.put(n, n.toString());
			}
			else if (n.toString().startsWith("DUD_"))
			{
				addDatasetCitation(n, "huang_benchmarking_2006",
						"http://pubs.acs.org/doi/abs/10.1021/jm0608356");
				addDatasetWeblink(n, "DUD-Directory", "http://dud.docking.org",
						"Directory of Useful Decoys");
				addDatasetCitation(n, "riniker_heterogeneous_2013",
						"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
				addDatasetWeblink(n, "Benchmarking-Platform",
						"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
				if (n == DataID.DUD_cdk2)
					datasetActivityDesc.put(n, "cyclin-dependent kinase");
				else if (n == DataID.DUD_hivrt)
					datasetActivityDesc.put(n, "HIV reverse transcriptase");
				else if (n == DataID.DUD_vegfr2)
					datasetActivityDesc.put(n, "vascular endothelial growth factor receptor");
				else
					throw new IllegalArgumentException();
				cfpDatasets.add(n);
				datasetCategory.put(n, VS_DATASETS);
				datasetSubCategory.put(n, "DUD");
			}
			else if (n.toString().startsWith("ChEMBL_"))
			{
				addDatasetCitation(n, "heikamp_large-scale_2011",
						"http://pubs.acs.org/doi/abs/10.1021/ci200199u");
				addDatasetCitation(n, "riniker_heterogeneous_2013",
						"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
				addDatasetWeblink(n, "Benchmarking-Platform",
						"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
				cfpDatasets.add(n);
				datasetCategory.put(n, VS_DATASETS);
				datasetSubCategory.put(n, "ChEMBL");
			}
			else if (n.toString().startsWith("MUV_"))
			{
				addDatasetCitation(n, "rohrer_maximum_2009",
						"http://pubs.acs.org/doi/abs/10.1021/ci8002649");
				addDatasetWeblink(n, "MUV-Datasets",
						"http://www.pharmchem.tu-bs.de/lehre/baumann/MUV.html",
						"Maximum Unbiased Validation (MUV) Datasets for Virtual Screening");
				addDatasetCitation(n, "riniker_heterogeneous_2013",
						"http://pubs.acs.org/doi/abs/10.1021/ci400466r");
				addDatasetWeblink(n, "Benchmarking-Platform",
						"https://github.com/rdkit/benchmarking_platform", "Benchmarking Platform");
				cfpDatasets.add(n);
				datasetCategory.put(n, VS_DATASETS);
				datasetSubCategory.put(n, "MUV");
			}
			//		n = "REID";
			//		datasetActivityDesc.put(n, "norm induction >= 11");
			//		sdfDatasets.put(n, "all_isoxazoles.sdf");
			//		sdfEndpoints.put(n, "activity");
			//		addDatasetCitation(n, "-", "-");
			//		addDatasetWeblink(n, "-", "-", "-");
			//		datasetCategory.put(n, OTHER_DATASETS);
			//		datasetSubCategory.put(n, n);
			//		for (int i : new int[] { 2, 3, 4, 11 })
			//		{
			//			n = "REID-" + i;
			//			datasetActivityDesc.put(n, "norm induction >= " + i);
			//			sdfDatasets.put(n, "reid-6.sdf");
			//			sdfEndpoints.put(n, "activity" + i);
			//			addDatasetCitation(n, "-", "-");
			//			addDatasetWeblink(n, "-", "-", "-");
			//			datasetCategory.put(n, OTHER_DATASETS);
			//			datasetSubCategory.put(n, n);
			//		}
			else if (n == DataID.LTKB)
			{
				datasetActivityDesc.put(n, "drug-induced liver injury (DILI)");
				addDatasetCitation(n, "-", "-");
				addDatasetWeblink(n, "Liver Toxicity Knowledge Base (LTKB)",
						"http://www.fda.gov/ScienceResearch/BioinformaticsTools/LiverToxicityKnowledgeBase/ucm2024036.htm",
						"Liver Toxicity Knowledge Base (LTKB)");
				datasetCategory.put(n, OTHER_DATASETS);
				datasetSubCategory.put(n, n.toString());
			}
		}

		//		datasetActivityDesc.put(Dataset.ChEMBL_101174,
		//				"pituitary adenylate cyclase-activating polypeptide type I receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_101395, "IgG receptor FcRn large subunit p51");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10102, "5-lipoxygenase activating protein");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10144, "bone morphogenetic protein 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12909, "ileal bile acid transporter");
		//		datasetActivityDesc.put(Dataset.ChEMBL_20130, "inhibitor of apoptosis protein");
		//		datasetActivityDesc.put(Dataset.ChEMBL_275, "retinoid X receptor alpha");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11061, "motilin receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10056, "DNA-dependent protein kinase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11096, "sodium/hydrogen exchanger 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10845, "phospholipase D1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11758, "glucagon-like peptide receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11402, "furin");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12725, "matriptase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_101219, "secreted frizzled-related protein 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_176, "Purinergic receptor P2Y12");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10087, "deoxycytidine kinase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_100098, "serine/threonine-protein kinase WEE1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10624, "serotonin 5a (5-HT5a) receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12659, "prostanoid DP receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10582, "cytosolic phospholipase A2");
		//		datasetActivityDesc.put(Dataset.ChEMBL_100862, "metastin receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_117, "somatostatin receptor 2");
		//		datasetActivityDesc.put(Dataset.ChEMBL_4,
		//				"voltage-gated T-type calcium channel alpha-1H subunit");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11635, "protein kinase C alpha");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11242, "Focal adhesion kinase 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_34, "fibronectin receptor beta");
		//		datasetActivityDesc.put(Dataset.ChEMBL_100077,
		//				"cell division cycle 7-related protein kinase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_193, "coagulation factor IX");
		//		datasetActivityDesc.put(Dataset.ChEMBL_80, "FK506-binding protein 1A");
		datasetActivityDesc.put(DataID.ChEMBL_165, "HERG");
		datasetActivityDesc.put(DataID.ChEMBL_10193, "carbonic anhydrase I");
		datasetActivityDesc.put(DataID.ChEMBL_15, "carbonic anhydrase II");
		datasetActivityDesc.put(DataID.ChEMBL_11489, "11-beta-hydroxysteroid dehydrogenase 1");
		datasetActivityDesc.put(DataID.ChEMBL_121, "serotonin transporter");
		datasetActivityDesc.put(DataID.ChEMBL_72, "dopamine D2 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_259, "cannabinoid CB2 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_10188, "MAP kinase p38 alpha");
		datasetActivityDesc.put(DataID.ChEMBL_108, "serotonin 2c (5-HT2c) receptor");
		datasetActivityDesc.put(DataID.ChEMBL_12952, "carbonic anhydrase IX");
		datasetActivityDesc.put(DataID.ChEMBL_93, "acetylcholinesterase");
		datasetActivityDesc.put(DataID.ChEMBL_10980,
				"vascular endothelial growth factor receptor 2");
		datasetActivityDesc.put(DataID.ChEMBL_19905, "melanin-concentrating hormone receptor 1");
		datasetActivityDesc.put(DataID.ChEMBL_107, "serotonin 2a (5-HT2a) receptor");
		datasetActivityDesc.put(DataID.ChEMBL_87, "cannabinoid CB1 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_17045, "cytochrome P450 3A4");
		datasetActivityDesc.put(DataID.ChEMBL_11140, "dipeptidyl peptidase IV");
		datasetActivityDesc.put(DataID.ChEMBL_114, "adenosine A1 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_90, "dopamine D4 receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_100166, "kinesin-like protein 1");
		datasetActivityDesc.put(DataID.ChEMBL_13001, "matrix metalloproteinase-2");
		datasetActivityDesc.put(DataID.ChEMBL_104, "monoamine oxidase B");
		datasetActivityDesc.put(DataID.ChEMBL_65, "cytochrome P450 19A1");
		datasetActivityDesc.put(DataID.ChEMBL_61, "muscarinic acetylcholine receptor M1");
		datasetActivityDesc.put(DataID.ChEMBL_10280, "histamine H3 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_51, "serotonin 1a (5-HT1a) receptor");
		datasetActivityDesc.put(DataID.ChEMBL_100, "norepinephrine transporter");
		datasetActivityDesc.put(DataID.ChEMBL_10260, "vanilloid receptor");
		datasetActivityDesc.put(DataID.ChEMBL_52, "alpha-2a adrenergic receptor");
		datasetActivityDesc.put(DataID.ChEMBL_11365, "cytochrome P450 2D6");
		datasetActivityDesc.put(DataID.ChEMBL_11359, "phosphodiesterase 4D");
		//		datasetActivityDesc.put(Dataset.ChEMBL_28, "thymidylate synthase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11536, "ghrelin receptor");
		datasetActivityDesc.put(DataID.ChEMBL_8, "tyrosine-protein kinase ABL");
		datasetActivityDesc.put(DataID.ChEMBL_10434, "tyrosine-protein kinase SRC");
		datasetActivityDesc.put(DataID.ChEMBL_12670, "tyrosine-protein kinase receptor FLT3");
		//		datasetActivityDesc.put(Dataset.ChEMBL_20014, "serine/threonine-protein kinase Aurora-A");
		//		datasetActivityDesc.put(Dataset.ChEMBL_234, "insulin-like growth factor I receptor");
		datasetActivityDesc.put(DataID.ChEMBL_12261, "c-Jun N-terminal kinase 1");
		datasetActivityDesc.put(DataID.ChEMBL_12209, "carbonic anhydrase XII");
		datasetActivityDesc.put(DataID.ChEMBL_25, "glucocorticoid receptor");
		datasetActivityDesc.put(DataID.ChEMBL_36, "progesterone receptor");
		datasetActivityDesc.put(DataID.ChEMBL_43, "beta-2 adrenergic receptor");
		datasetActivityDesc.put(DataID.ChEMBL_219, "muscarinic acetylcholine receptor M3");
		datasetActivityDesc.put(DataID.ChEMBL_130, "dopamine D3 receptor");
		datasetActivityDesc.put(DataID.ChEMBL_105, "serotonin 1d (5-HT1d) receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11336, "neuropeptide Y receptor type 5");
		//		datasetActivityDesc.put(Dataset.ChEMBL_20174, "G protein-coupled receptor");
		datasetActivityDesc.put(DataID.ChEMBL_126, "cyclooxygenase-2");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11225, "renin");
		datasetActivityDesc.put(DataID.ChEMBL_12252, "beta-secretase 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11682, "glycine transporter 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_134, "vasopressin V1a receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_116, "oxytocin receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11265, "somatostatin receptor 5");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10475, "neuropeptide Y receptor type 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12679, "C5a anaphylatoxin chemotactic receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10579, "C-C chemokine receptor type 4");
		datasetActivityDesc.put(DataID.ChEMBL_11575, "C-C chemokine receptor type 2");
		//		datasetActivityDesc.put(Dataset.ChEMBL_18061,
		//				"sodium channel protein type IX alpha subunit");
		//		datasetActivityDesc.put(Dataset.ChEMBL_237, "leukotriene A4 hydrolase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_276, "phosphodiesterase 4A");
		datasetActivityDesc.put(DataID.ChEMBL_11534, "cathepsin S");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10198,
		//				"voltage-gated potassium channel subunit Kv1.5");
		datasetActivityDesc.put(DataID.ChEMBL_10498, "cathepsin L");
		datasetActivityDesc.put(DataID.ChEMBL_12911, "cytochrome P450 2C9");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12968, "orexin receptor 2");
		datasetActivityDesc.put(DataID.ChEMBL_100579, "nicotinic acid receptor 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_100126, "serine/threonine-protein kinase B-raf");
		datasetActivityDesc.put(DataID.ChEMBL_10378, "cathepsin B");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10417, "P2X purinoceptor 7");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10752,
		//				"inhibitor of nuclear factor kappa B kinase beta subunit");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10773, "interleukin-8 receptor B");
		datasetActivityDesc.put(DataID.ChEMBL_11631, "sphingosine 1-phosphate receptor Edg-1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_10927, "urotensin II receptor");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11085, "melatonin receptor 1B");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11442, "liver glycogen phosphorylase");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11279, "metabotropic glutamate receptor 1");
		//		datasetActivityDesc.put(Dataset.ChEMBL_11488, "estradiol 17-beta-dehydrogenase 3");
		//		datasetActivityDesc.put(Dataset.ChEMBL_12840,
		//				"macrophage colony stimulating factor receptor");		

		datasetActivityDesc.put(DataID.MUV_466, "S1P1 rec. (GPCR) Agonist");
		datasetActivityDesc.put(DataID.MUV_548, "PKA (Kinase) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_600, "SF1 (Nuclear Receptor) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_644, "Rho-Kinase2 Inhibitor");
		datasetActivityDesc.put(DataID.MUV_652, "HIV RT-RNase Inhibitor");
		datasetActivityDesc.put(DataID.MUV_689, "Eph rec. A4 (Rec. Tyr. Kinase) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_692, "SF1 (Nuclear Receptor) Agonist");
		datasetActivityDesc.put(DataID.MUV_712, "HSP 90 (Chaperone) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_713, "ER-a-Coact. Bind. (PPI) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_733, "ER-ß-Coact. Bind. (PPI) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_737, "ER-a-Coact. Bind. (PPI) Potentiator");
		datasetActivityDesc.put(DataID.MUV_810, "FAK (Kinase) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_832, "Cathepsin G (Protease) Inhibitor");
		//		datasetActivityDesc.put(Dataset.MUV_846, "FXIa (Protease) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_852, "FXIIa (Protease) Inhibitor");
		datasetActivityDesc.put(DataID.MUV_858, "D1 rec. (GPCR) Allosteric Modulator");
		datasetActivityDesc.put(DataID.MUV_859, "M1 rec. (GPCR) Allosteric Modulator");
	}

	private static KeyValueFileStore<String, CDKDataset> parsedDatasets;

	public static CDKDataset getDataset(DataID dataset)
	{
		if (!datasets.containsKey(dataset))
		{
			try
			{
				if (parsedDatasets == null)
					parsedDatasets = new KeyValueFileStore<>(DATA_FOLDER + "/parsed", false, true,
							null, true);
				CDKDataset data;
				if (!parsedDatasets.contains(dataset.toString()))
				{
					if (sdfDatasets.containsKey(dataset))
						data = getDatasetFromSDF(dataset);
					else if (new File(DATA_FOLDER + File.separator + dataset + ".csv").exists())
						data = getDatasetFromCSV(dataset);
					else if (new File(DATA_FOLDER + File.separator + dataset + ".smi").exists())
						data = getDatasetFromSMILES(dataset);
					else
						throw new RuntimeException(
								"dataset not found: " + dataset + " in " + DATA_FOLDER);
					parsedDatasets.store(dataset.toString(), data);
				}
				else
					data = parsedDatasets.get(dataset.toString());
				datasets.put(dataset, data);
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}
		}
		return datasets.get(dataset);
	}

	private static CDKDataset getDatasetFromSDF(DataID dataset)
			throws FileNotFoundException, CDKException, IOException
	{
		return DataLoader.getDatasetFromSDF(dataset.toString(),
				DATA_FOLDER + File.separator + sdfDatasets.get(dataset), sdfEndpoints.get(dataset));
	}

	private static CDKDataset getDatasetFromCSV(DataID dataset) throws CDKException
	{
		return DataLoader.getDatasetFromCSV(dataset.toString(),
				DATA_FOLDER + File.separator + dataset + ".csv", 0, 1);
	}

	private static CDKDataset getDatasetFromSMILES(DataID dataset) throws CDKException
	{
		return DataLoader.getDatasetFromSMILES(dataset.toString(),
				DATA_FOLDER + File.separator + dataset + ".smi");
	}

	public static boolean hasDuplicates(String name)
	{
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

	private static boolean hasDuplicatesSDF(String name) throws Exception
	{
		List<String> smiles = new ArrayList<>();
		ISimpleChemObjectReader reader = new ReaderFactory().createReader(new InputStreamReader(
				new FileInputStream(DATA_FOLDER + File.separator + sdfDatasets.get(name))));
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

	private static boolean hasDuplicatesCSV(String name) throws CDKException
	{
		List<String> smiles = new ArrayList<>();
		CSVFile csv = FileUtil.readCSV(DATA_FOLDER + File.separator + name + ".csv");
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

	public static Set<DataID> cfpDatasets()
	{
		return cfpDatasets;
	}

	public static DataID[] cfpDatasetsSorted()
	{
		return cfpDatasets.stream().sorted(CFPDataComparator).toArray(size -> new DataID[size]);
	}

	private static Set<DataID> cfpDatasetsFiltered(Predicate<DataID> filter)
	{
		return cfpDatasets.stream().filter(filter).collect(Collectors.toSet());
	}

	public static Set<DataID> cfpDatasetsCPDB()
	{
		return cfpDatasetsFiltered(d -> d.toString().startsWith("CPDBAS_"));
	}

	public static Set<DataID> cfpDatasetsChEMBL()
	{
		return cfpDatasetsFiltered(d -> d.toString().startsWith("ChEMBL_"));
	}

	public static Set<DataID> cfpDatasetsMUV()
	{
		return cfpDatasetsFiltered(d -> d.toString().startsWith("MUV_"));
	}

	public static Set<DataID> cfpDatasetsDUD()
	{
		return cfpDatasetsFiltered(d -> d.toString().startsWith("DUD_"));
	}

	public static Set<DataID> cfpDatasetsBalanced()
	{
		return cfpDatasetsCategory(BALANCED_DATASETS);
	}

	public static Set<DataID> cfpDatasetsCategory(String category)
	{
		return cfpDatasetsFiltered(d -> datasetCategory.get(d).equals(category));
	}

	public static Set<DataID> cfpDatasetsSubCategory(String subCategory)
	{
		return cfpDatasetsFiltered(d -> datasetSubCategory.get(d).equals(subCategory));
	}

	//	public boolean exists(String name)
	//	{
	//		if (datasets.containsKey(name))
	//			return true;
	//		if (sdfDatasets.containsKey(name))
	//			return new File(dataFolder + File.separator + sdfDatasets.get(name)).exists();
	//		return new File(dataFolder + File.separator + name + ".csv").exists();
	//	}

	public static ResultSet getInfo(boolean cite, DataID... datasets)
	{
		ResultSet set = new ResultSet();
		set.setNicePropery("size", "compounds");
		for (DataID dataset : datasets)
		{
			//			System.out.println(n);

			int rIdx = set.addResult();
			set.setResultValue(rIdx, "category", datasetCategory.get(dataset));
			set.setResultValue(rIdx, "name", dataset.toString());//.replaceAll("_", " "));
			CDKDataset d = getDataset(dataset);
			String classV[] = DataLoader.getClassValues(d.endpoints);
			int activeIdx = DataLoader.getActiveIdx(classV);
			set.setResultValue(rIdx, "size", d.getSmiles().size());
			CountedSet<String> endp = CountedSet.create(d.getEndpoints());
			set.setResultValue(rIdx, "active", endp.getCount(classV[activeIdx]));
			set.setResultValue(rIdx, "in-active", endp.getCount(classV[1 - activeIdx]));
			//			set.setResultValue(rIdx, "#in-active", endp.getCount(classV[1 - activeIdx]));
			//			set.setResultValue(rIdx, "activity", ArrayUtil.toString(classV));
			set.setResultValue(rIdx, "target", datasetActivityDesc.get(dataset));

			//			set.setResultValue(rIdx, "dataset-weblink", CollectionUtil.toString(datasetWeblinks.get(n)));

			String source;
			if (cite)
			{
				String cit = "\\cite{";
				for (Source s : citation.get(dataset))
					if (!(s instanceof WebSource))
						cit += s.citationKey + ",";
				cit = cit.substring(0, cit.length() - 1) + "}";
				source = cit;
			}
			else
			{
				source = "";
				for (Source s : citation.get(dataset))
					if (s instanceof WebSource)
						source += "\\url{" + s.citationUrl + "}, ";
				source = source.substring(0, source.length() - 2);
			}
			set.setResultValue(rIdx, "source", source);
		}
		return set;
	}

	public static ResultSet getCategoryInfo(boolean cite, DataID... dataset)
	{
		ResultSet set = getInfo(cite, dataset);

		for (int idx = 0; idx < set.getNumResults(); idx++)
		{
			DataID n = DataID.valueOf(set.getResultValue(idx, "name").toString());
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

	public static String getDatasetEndpoint(DataID d)
	{
		return datasetActivityDesc.get(d);
	}

	public static Set<String> getDatasetURLs(DataID d)
	{
		HashSet<String> urls = new LinkedHashSet<>();
		for (Source s : citation.get(d))
			urls.add(s.citationUrl);
		return urls;
	}

	public static Map<String, String> getModelDatasetCitations(DataID d)
	{
		HashMap<String, String> map = new LinkedHashMap<>();
		for (Source s : citation.get(d))
			map.put(s.getCitationName(), s.citationUrl);
		return map;
	}

	public boolean isVirtualScreeningDataset(String name)
	{
		return datasetCategory.get(name) == VS_DATASETS;
	}

	public static void main(String[] args)
	{
		{
			System.out.println(getCategoryInfo(true, cfpDatasetsSorted()).toNiceString());
		}
		{
			DataID data[] = DataID.values();
			//		Arrays.sort(data, CFPDataComparator);
			//		System.out.println(ArrayUtil.toString(data));
			data = new DataID[] { DataID.CPDBAS_Mouse };
			for (DataID dat : data)
			{
				CDKDataset da = DataProvider.getDataset(dat);
				if (!da.warnings.isEmpty())
				{
					System.out.println(dat + " Warnings:");
					for (String warn : da.warnings)
						System.out.println("* " + warn);
					System.out.println();
				}
			}
		}
		{
			System.out.println(getInfo(true, cfpDatasetsSorted()).toNiceString());
		}
	}

}
