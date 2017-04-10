import java.io.IOException;
import java.io.*;
import java.nio.file.*;
import java.util.*;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.*;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;


// To get the "genoCode" from a bam file: change the fillHelper function

public class GenoTable 
{
	private String d_genoVCF;
	private String d_SNPfile;
    private String d_cellNames;
    private String d_chr = "1";
    
    private int d_numbSamples;
    private int d_numbCells;
    private boolean d_loopOverAllChr;
    
    private RandomAccessGenotypeData d_dnaGenotypes; 
    private Hashtable<String, Scores> d_genoTable;
    private double[] d_alphas;
    private Hashtable<String, double[]> d_doublets;
    private int d_tVal; // Cut-off value for (doublet - singlet) likelihoods
	
    public GenoTable(String genoVCF, String SNPfile, String cellNameFile, boolean allChr, int alpha, int tVal)
    {
    	d_genoVCF = genoVCF;
    	d_SNPfile = SNPfile;
    	d_cellNames = cellNameFile;
    	d_loopOverAllChr = allChr;
    	d_tVal = tVal;
    	parseAlpha(alpha);
    }
    
    private void parseAlpha(int alpha)
    {
    	if (alpha == 1)
    		d_alphas = new double[] {0.5};
    	else if (alpha == 3)
    		d_alphas = new double[] {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
    	else if (alpha == 4)
    		d_alphas = new double[] {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
    	else 
    		d_alphas = new double[] {0.25, 0.5, 0.75};	
    }
    
    public void run() throws IOException
    {
    	createTable();
    	fillTable();
    
    }
    
    public void createTable() throws IOException
    {
    	// Read DNA Genotypes
		try 
		{
			d_dnaGenotypes = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(d_genoVCF);
			
		} 
		catch (IOException ex) 
		{
			System.err.println("Error reading genotype VCF file (from the -g option).\n");
			throw (ex);
		}
		
		d_numbSamples = d_dnaGenotypes.getSampleNames().length;
		System.out.println("The number of samples in the DNA genotypes is: " + d_numbSamples);
		
		// Make hashtable, for each cell a list of d_numbSamples scores
		d_genoTable = new Hashtable<String, Scores>();
		d_doublets = new Hashtable<String, double[]>();
		
		// Read cell names
		try
		{
			for (String fline : Files.readAllLines(Paths.get(d_cellNames))) 
			{    
				d_genoTable.put(fline, new Scores(d_numbSamples));
				d_doublets.put(fline, new double[(d_numbSamples * (d_numbSamples - 1) / 2) * (d_alphas.length)]);

			}
			d_numbCells = d_genoTable.size();
			System.out.println("The number of cell barcodes is: " + d_numbCells + "\n");
		}
		catch(IOException ex) 
		{
			System.err.println("Error reading cell name file (from the -n option).\n");
			throw (ex);
		}
    }
    
    public void fillTable() throws IOException
    {
    	System.out.println("Started assigning SNPs per cell to samples...");
    	
    	if (d_loopOverAllChr)
    	{	
    		System.out.println("Looping over all chromosomes (1 to 22).\n");
    		for (int chr = 1; chr < 23; ++chr)
    		{
    			System.out.println("Starting chromosome " + chr);    			
    			d_SNPfile = replaceChr(d_SNPfile, chr);
    			d_genoVCF = replaceChr(d_genoVCF, chr);    			
    			d_chr = String.valueOf(chr);  // must be placed AFTER replaceChr() calls!
	
    			// Reload d_dnaGenotypes for the new CHR
    			String[] oldSampleNames = d_dnaGenotypes.getSampleNames();
    			try 
    			{
    				d_dnaGenotypes = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(d_genoVCF);
    			} 
    			catch (IOException ex) 
    			{
    				System.err.println("Error reading genotype VCF file (from the -g option).\n");
    				throw (ex);
    			}
    			
    			// Check if sample names are the same between current and previous chromosome of genoVCF
    			String[] newSampleNames= d_dnaGenotypes.getSampleNames();
    			boolean sameNames = true;
    			for (int idx = 0; idx < d_numbSamples; ++idx)
    			{
    				if(!Arrays.asList(newSampleNames).contains(oldSampleNames[idx]))
    					sameNames = false;
    			}
    			
    			if (!(newSampleNames.length == d_numbSamples) || !sameNames)
    			{
    				System.out.println("WARNING! Different sample names or number of samples in the genoVCF of chr " + chr + " and the previous chr. Not continuing.");
    				return;
    			}

    			fillHelper();
    		}
    	}
    	else
    	{ 
    		fillHelper();
    	}
    } 
    
    private String replaceChr(String str, int chr)
    {
    	String oldChr = "chr" + d_chr;
    	String newChr = "chr" + chr;
    	String newStr = str.replaceAll(oldChr, newChr);
    	return newStr;
    }
    
    private void fillHelper() throws IOException
    {
       	int pos;
    	int[] genoCode = new int[4]; // in alphabetic order: A, C, G, T
    	String cell_id;
    	
    	// parse line from d_SNPfile, get pos, cell_id and reads for the alleles
    	BufferedReader reader = null;
    	try {
    		String curLine;
    		reader = new BufferedReader(new FileReader(d_SNPfile));
    		int processedlines = 0;
    		while ((curLine = reader.readLine()) != null)
    		{
    			String[] strArr = curLine.split("\t");
    			if (!strArr[0].equals("CHROM"))
    			{	
    				d_chr = strArr[0];
	    			cell_id = strArr[4];
	    			//System.out.println(cell_id);
	    			pos = Integer.parseInt(strArr[1]);
	    			for (int idx = 0; idx < 4; ++idx)
	    				genoCode[idx] = Integer.parseInt(strArr[5 + idx]);

	    		   	// Match to vcf	
	    		  	match(cell_id, pos, genoCode);
	    		  	++processedlines;
	    		  	if (processedlines % 25000 == 0)
	    		  		System.out.println(processedlines + " lines processed...");
    			}
    			
 
    		}
    	}
    	catch (IOException e) {
    		System.err.println("Error reading SNP file (from the -s option).\n");
    		e.printStackTrace();
        } 
    	finally 
    	{
            try 
            {
                if (reader != null)
                	reader.close();
            } 
            catch (IOException ex) 
            {
                ex.printStackTrace();
            }
        }

    	System.out.println("Done!\n");
    }
    
    
    private void match(String cell_id, int pos, int[] genoCode)
    {
    	// check if cell_id matches
    	if (!d_genoTable.containsKey(cell_id))
    		return;

    	GeneticVariant snp = d_dnaGenotypes.getSnpVariantByPos(d_chr,pos);
    	
    	if (!(snp != null))
    		return; 
    	
    	if (!snp.isBiallelic())
    		return;
    	
    	// Check that all samples have the GP field filled?
    	Allele ref = snp.getRefAllele();
    	Allele alt = snp.getAlternativeAlleles().get(0);
		float[][] snpGP = snp.getSampleGenotypeProbilities();

		// Check that all samples have the GP field filled
		boolean allFilled = true;
		for (int samp = 0; samp < d_numbSamples; ++samp)
		{
			if (snpGP[samp][0] + snpGP[samp][1] + snpGP[samp][2] < 1.0)
				allFilled = false;
		}
		
		if (!allFilled)
			return;
		
		// If you have to convert PL to GP? Need MAF?
		//double maf = snp.getMinorAlleleFrequency();
		
		// Do comparison (in class SNPCompare)
		SNPCompare comp = new SNPCompare(genoCode, snpGP, ref, alt, d_numbSamples, d_alphas);
		double[] result = comp.run();
		
		if (!comp.goodSNP())
			return;
		
		Scores score = d_genoTable.get(cell_id);
		
		score.addScores(result);
		score.increaseNSnps();
		
		d_genoTable.put(cell_id, score);
		
		
		double[] doubletResult = comp.runDoublets();
		double[] doubletOldArr = d_doublets.get(cell_id);
		
		for (int idx = 0; idx < (d_numbSamples * (d_numbSamples - 1) / 2) * (d_alphas.length); ++idx)
			doubletOldArr[idx] += doubletResult[idx];
		
		d_doublets.put(cell_id, doubletOldArr);	
    }
   
    
    public int getSamp1(int combi)
    {
        ++combi; 
    	int combiNew = 0;
    	for (int idx = 0; idx < d_numbSamples; ++idx )
    		for (int idx2 = idx + 1; idx2 < d_numbSamples; ++idx2)
    		{
    			++combiNew;
    			if (combiNew == combi)
    				return idx;
    		}
    	return -1;
    }
    
    public int getSamp2(int combi)
    {
    	++combi;
    	int combiNew = 0;
    	for (int idx = 0; idx < d_numbSamples; ++idx )
    		for (int idx2 = idx + 1; idx2 < d_numbSamples; ++idx2)
    		{
    			++combiNew;
    			if (combiNew == combi)
    				return idx2;
    		}
    	return -1;
    }
    
    public void writeTable(String outputPath) throws IOException
    {
    	System.out.println("Writing table to: " + outputPath);
    	
    	// Make file
    	// Don't use ~ in outputPath
    	try 
    	{

  	      File file = new File(outputPath);

  	      if (file.createNewFile())
  	        System.out.println("File is created.");
  	      else
  	        System.out.println("File already exists. Overwritten.");

      	} 
    	catch (IOException e) 
    	{
  	      e.printStackTrace();
      	}
    	
    	// Write to file
    	try{
    	    PrintWriter writer = new PrintWriter(outputPath, "UTF-8");

    	    //writer.println("\t" + Arrays.toString(d_dnaGenotypes.getSampleNames()).replaceAll("[\\[ | \\]]", "").replaceAll(",", "\t") + "\t" + "nSNPs_tested");
    	    writer.println("cell_id\tSinglet_samp\tDoub_samp1\tDoub_samp2\talpha\tllkDoublet-llkSinglet\tllkDoublet\tllkSingletSamp1\tllkSingletSamp2\tnSNPs_tested\toutcome\tassigned_sample(s)");
    	    
    	    for (String cell_id : Files.readAllLines(Paths.get(d_cellNames))) 
    		{
    	    	//int maxAt = maxValueAt(d_genoTable.get(fline).getScores());
    	    	//writer.println(fline + "\t" + Arrays.toString(d_genoTable.get(fline).getScores()).replaceAll("[\\[ | \\]]", "").replaceAll(",", "\t") + "\t" + d_genoTable.get(fline).getNSnps() + "\t" + d_dnaGenotypes.getSampleNames()[maxAt]);
    	    	//writer.println(fline + "\t" + Arrays.toString(d_doublets.get(fline)).replaceAll("[\\[ | \\]]", "").replaceAll(",", "\t") + "\t");
    	    	int maxDoubAt = maxValueAt(d_doublets.get(cell_id));
    	    	
    	    	int combi = maxDoubAt / (d_alphas.length);
    	    	int alpha = maxDoubAt - (d_alphas.length) * (maxDoubAt / (d_alphas.length));
    	    	int maxCellAt = maxValueAt(d_genoTable.get(cell_id).getScores());
    	    	    	
    	    	double maxDoub = d_doublets.get(cell_id)[maxDoubAt];
    	    	double maxCell = d_genoTable.get(cell_id).getScores()[maxCellAt];
    	    	double max = maxDoub - maxCell;

    	    	if (max > d_tVal)
    	     		writer.println(cell_id + "\t" + maxCellAt + "\t" + getSamp1(combi) + "\t" + getSamp2(combi)+ "\t" + d_alphas[alpha] + "\t" + 
    	     				max + "\t" + maxDoub + "\t" + d_genoTable.get(cell_id).getScores()[getSamp1(combi)] + "\t" + 
    	     				d_genoTable.get(cell_id).getScores()[getSamp2(combi)] + "\t" + d_genoTable.get(cell_id).getNSnps() +"\tDoublet\t" + d_dnaGenotypes.getSampleNames()[getSamp1(combi)] + 
    	     				"," + d_dnaGenotypes.getSampleNames()[getSamp2(combi)]);
    	    	else if (max < -d_tVal)
    	    		writer.println(cell_id + "\t" + maxCellAt + "\t" + getSamp1(combi) + "\t" + getSamp2(combi)+ "\t" + d_alphas[alpha] + "\t" + 
    	     				max + "\t" + maxDoub + "\t" + d_genoTable.get(cell_id).getScores()[getSamp1(combi)] + "\t" + 
    	     				d_genoTable.get(cell_id).getScores()[getSamp2(combi)] + "\t" + d_genoTable.get(cell_id).getNSnps() + "\tSinglet\t" + d_dnaGenotypes.getSampleNames()[maxCellAt]);
    	    	else
    	    		writer.println(cell_id + "\t" + maxCellAt + "\t" + getSamp1(combi) + "\t" + getSamp2(combi)+ "\t" + d_alphas[alpha] + "\t" + 
    	     				max + "\t" + maxDoub + "\t" + d_genoTable.get(cell_id).getScores()[getSamp1(combi)] + "\t" + 
    	     				d_genoTable.get(cell_id).getScores()[getSamp2(combi)] + "\t" + d_genoTable.get(cell_id).getNSnps() +"\tInconclusive");
   	
    		}
    	    writer.close();
    	} 
    	catch (IOException e) 
    	{
    	   System.out.println("No table could be written. Error message:" + e.getMessage()) ;
    	}	
    	System.out.println("Done!\n");
    }
    
    
    private static int maxValueAt(double[] doubs) 
    {
        // Assuming a unique maximum...
    	double max = doubs[0];
        int maxAt = 0;
        for (int ktr = 0; ktr < doubs.length; ++ktr) {
            if (doubs[ktr] > max) {
                max = doubs[ktr];
                maxAt = ktr;
            }
        }
        return maxAt;
    }
     
 }
