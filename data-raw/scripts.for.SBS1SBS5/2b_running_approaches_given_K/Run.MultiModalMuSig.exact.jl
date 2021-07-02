
# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH <- string("<path_to_results_on_SBS1-SBS5-correlated_datasets>")
#
# cd(PATH)

using MultiModalMuSig
using CSV
using DataFrames
using VegaLite
using Random

## Specify top level directories of spectra dataset and final output
topLevelFolder4Data = "../research_data/0.Input_datasets"
topLevelFolder4Run = "../research_data/2b.Full_output_K_as_2"



## Specify datasetnames
slopes = (0.1,0.5,1,2,10);
Rsqs = (0.1,0.2,0.3,0.6);
datasetNames = [];

for slope in slopes
    for Rsq in Rsqs
        push!(datasetNames,string("S.",slope,".Rsq.",Rsq))
    end
end

## Specify seedsInUse
seedsInUse = [1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753];



## Cycle for each dataset to run MultiModalMuSig.LDA
for seedInUse in seedsInUse
    for datasetName in datasetNames
        outDir = paste0(topLevelFolder4Run,"/MultiModalMuSig.LDA.results/",datasetName,"/seed.",seedInUse)
        snv_counts = CSV.read(outDir*"/ground.truth.syn.catalog.tsv", delim='\t');

        ## Run model fitting
        ## After the model is fitted, we may obtain extracted signatures and inferred exposures
        X = format_counts_lda(snv_counts);

        ## Specify the best K
        KBest = 2

        ## Obtain extracted signatures given the best K.
        ## In MultiModalMuSig format, later converted to ICAMS format
        Random.seed!(seedInUse);
        model = LDA(KBest, 0.1, 0.1, 96, X);
        fit!(model, tol=1e-5);
        snv_signatures = DataFrame(model.β);
        sig_names = names(snv_signatures);
        snv_signatures[:term] = snv_counts[:term];
        snv_signatures = snv_signatures[:,[:term;sig_names]]; # Shuffle the columns of signature DataFrame

        # You can't change row names by subsetting!
        CSV.write(outDir*"/extracted.signatures.tsv",
        snv_signatures,
        delim = "\t");

        ## Obtain inferred exposures
        snv_exposures = DataFrame(model.θ);
        sample_names = names(snv_counts)[2:(size(snv_counts)[2])];
        names!(snv_exposures, sample_names);# Set column names for DataFrame snv_exposures
        # For each tumor, Multiply exposure proportion with the total number of mutations
        for sample_name in sample_names
            snv_exposures[:,sample_name] = snv_exposures[:,sample_name] * sum(snv_counts[:,sample_name]);
        end
        # Add signature names.
        # Julia DataFrame does not support row names to be "character"!
        snv_exposures[:signature_name] = sig_names
        snv_exposures = snv_exposures[:,[:signature_name; sample_names]]

        # Write inferred exposures
        CSV.write(outDir*"/inferred.exposures.tsv",
        snv_exposures,
        delim = "\t")
    end
end













## Cycle for each dataset to run MultiModalMuSig.CTM
for seedInUse in seedsInUse
    for datasetName in datasetNames
        outDir = paste0(topLevelFolder4Run,"/MultiModalMuSig.CTM.results/",datasetName,"/seed.",seedInUse)
        snv_counts = CSV.read(outDir*"/ground.truth.syn.catalog.tsv", delim='\t');

        ## Run model fitting
        ## After the model is fitted, we may obtain extracted signatures and inferred exposures
        X = format_counts_mmctm(snv_counts);

        ## Specify the best K
        KBest = [2];

        ## Obtain extracted signatures given the best K.
        ## In MultiModalMuSig format, later converted to ICAMS format
        Random.seed!(seedInUse);
        model = MMCTM(KBest, [0.1], X); ## You cannot write [KBest] here!
        ## as KBest is already an array, [KBest] would be an nested array!
        fit!(model, tol=1e-5);
        snv_signatures = DataFrame(hcat(model.ϕ[1]...));
        sig_names = names(snv_signatures);
        snv_signatures[:term] = snv_counts[:term];
        snv_signatures = snv_signatures[:,[:term;sig_names]]; # Shuffle the columns of signature DataFrame

        # You can't change row names by subsetting!
        CSV.write(outDir*"/extracted.signatures.tsv",
        snv_signatures,
        delim = "\t");

        ## Obtain inferred exposures
        snv_exposures = DataFrame(vcat(model.props...));
        sample_names = names(snv_counts)[2:(size(snv_counts)[2])];
        names!(snv_exposures, sample_names);# Set column names for DataFrame snv_exposures
        # For each tumor, Multiply exposure proportion with the total number of mutations
        for sample_name in sample_names
            snv_exposures[:,sample_name] = snv_exposures[:,sample_name] * sum(snv_counts[:,sample_name]);
        end
        # Add signature names.
        # Julia DataFrame does not support row names to be "character"!
        snv_exposures[:signature_name] = sig_names
        snv_exposures = snv_exposures[:,[:signature_name; sample_names]]

        # Write inferred exposures
        CSV.write(outDir*"/inferred.exposures.tsv",
        snv_exposures,
        delim = "\t")
    end
end
