using MultiModalMuSig
using CSV
using DataFrames
using VegaLite
using Random


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

## Set Eta (prior for signature shape beta) to 0.1 to 1.
## setting below 1 may prefer spiky signatures, 
## and thus may cause oversplit.
## set alpha (prior for exposure theta) to 0.1.
alpha = 0.1;
etas = [2.5, 5.0, 7.5, 10.0,
            15.0, 20.0, 30.0, 50.0, 75.0];



## Cycle for each dataset to run MultiModalMuSig.LDA
for eta in etas
    for seedInUse in seedsInUse
        for datasetName in datasetNames

            snv_counts = CSV.read(datasetName*"/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha."*string(alpha)*".eta."*string(eta)*".results/seed."*
            string(seedInUse)*"/ground.truth.syn.catalog.tsv", delim='\t');

            ## Run model fitting
            ## After the model is fitted, we may obtain extracted signatures and attributed exposures
            X = format_counts_lda(snv_counts);

            ## Save the likelihood for each K.
            likelihoods = [];
            for K in 2:10
                ## Specify seed used before fitting model for each dataset,
                ## for each signature number (K)
                Random.seed!(seedInUse);
                model = LDA(K, alpha, eta, 96, X);
                fit!(model, tol=1e-5);
                likelihood = MultiModalMuSig.calculate_loglikelihood(model);
                append!(likelihoods, likelihood);
            end

            ## Find which K is the best likelihood.
            max_likelihood = maximum(likelihoods);
            KBest = [i for (i, x) in enumerate(likelihoods) if x == max_likelihood];
            KBest = KBest[1]; ## LDA accepts K as an integer rather than an array

            ## Obtain extracted signatures given the best K.
            ## In MultiModalMuSig format, later converted to ICAMS format
            Random.seed!(seedInUse);
            model = LDA(KBest, alpha, eta, 96, X);
            fit!(model, tol=1e-5);
            snv_signatures = DataFrame(model.β);
            sig_names = names(snv_signatures);
            snv_signatures[:term] = snv_counts[:term];
            snv_signatures = snv_signatures[:,[:term;sig_names]]; # Shuffle the columns of signature DataFrame

            # You can't change row names by subsetting!
            CSV.write(datasetName*"/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha."*string(alpha)*".eta."*string(eta)*".results/seed."*
            string(seedInUse)*"/extracted.signatures.tsv",
            snv_signatures,
            delim = "\t");

            ## Obtain attributed exposures
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

            # Write attributed exposures
            CSV.write(datasetName*"/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha."*string(alpha)*".eta."*string(eta)*".results/seed."*
            string(seedInUse)*"/inferred.exposures.tsv",
            snv_exposures,
            delim = "\t")
        end
    end
end










