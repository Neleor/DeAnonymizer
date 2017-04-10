

import java.io.IOException;
import org.apache.commons.cli.*;


public class Main 
{
	public static void main(String[] args) throws IOException 
	{
	    Options options = new Options();
	    options.addOption("g", "genotypeVcf", true, "Input genotype file path");
	    options.addOption("s", "SNPfile", true, "Input SNP-reads file from countUMI-script, of matching chromosome of genotypeVCF");
	    options.addOption("n", "cellnames", true, "File containing names of all cells to be considered");
	    options.addOption("o", "output", true, "Output path + file");
	    options.addOption("a", "alpha", true, "OPTIONAL. Value (1-4), where 1: {0.5}, 2: {0.25, 0.5, 0.75}, 3: {0.10, 0.20, ..., 0.90}, 4: {0.05, 0.10, ..., 0.95}, default = 2.");
	    options.addOption("P", "parseAllChr", false, "OPTIONAL. When added, the genotypeVCF and SNPfile names are assumed to have \"chr1\" in them uniquely, where can be looped over, up to chr22. Assumed to start at chr 1.");
	    options.addOption("t", "tValCutoff", false, "OPTIONAL. Determine tValCutoff (for doublet llk - singlet llk). Default = 1.");
		   
	    // use DefaultParser(); when using apache commons > 1.3, else use GnuParser(); (Deprecated)
	    CommandLineParser parser = new DefaultParser();
	    
	    try 
	    {
            CommandLine line = parser.parse(options, args);
            if(!line.hasOption("g") || !line.hasOption("s") || !line.hasOption("o")) 
            {
                HelpFormatter formatter = new HelpFormatter();
                String header = "DeAnonymizer. \nThe program will compute the most likely sample to which each input cell belongs."
                		+ " Provide as input: genotype VCF of the samples to be considered + list of cell-ids + a file with (SNP position, cell-id, genotype), "
                		+ "as generated with 1_count_umi.py.\n"
                		+ "Enter the output path WITHOUT a ~ in the path.\n\n";
                formatter.printHelp(header, options, true);
                System.exit(1);
            }
            
            String genoVCF = line.getOptionValue("g");
            String SNPfile = line.getOptionValue("s");
            String cellNameFile = line.getOptionValue("n");
            String outputPath = line.getOptionValue("o");
            
            int alpha = 2;
            if(line.hasOption("a"))
            	alpha = Integer.parseInt(line.getOptionValue("a"));
            
            int tVal = 1;
            if(line.hasOption("t"))
            	alpha = Integer.parseInt(line.getOptionValue("t"));
            
            GenoTable genoTable;
                       
            // Iterate over all (22) chr or not?
	        if (line.hasOption("P"))
	        {
	          	genoTable = new GenoTable(genoVCF, SNPfile, cellNameFile, true, alpha, tVal);
	          	genoTable.run();
	        }
	        else 
	        {
	            genoTable = new GenoTable(genoVCF, SNPfile, cellNameFile, false, alpha, tVal);
	            genoTable.run(); 
	        }
   
            genoTable.writeTable(outputPath);

            
            System.out.println("Finished the program!");
            

	    }
	    catch (ParseException exp) 
	    {
            System.out.println("Parse exception:" + exp.getMessage());
        } 
	}
}
	    


	





