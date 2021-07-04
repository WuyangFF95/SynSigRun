
# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH = string("<path_to_results_on_SBS1-SBS5-correlated_datasets")
#
# cd(PATH)

using MultiModalMuSig
using CSV
using DataFrames
using NamedArrays
using VegaLite
using Random

## Specify top level directories of spectra dataset and final output
topLevelFolder4Data = "./0.Input_datasets"
topLevelFolder4Run = "./2a.Full_output_K_unspecified"



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
seedsInUse = seedsInUse[1:5];


## Cycle for each dataset to run MultiModalMuSig.LDA
for seedInUse in seedsInUse
    for datasetName in datasetNames
        outDir = string(topLevelFolder4Run,"/MultiModalMuSig.LDA.results/",datasetName,"/seed.",seedInUse)
        snv_counts = CSV.read(outDir*"/ground.truth.syn.catalog.tsv", delim='\t');
        
        print("\n\n==================================\n")
        print(string("Running MultiModalMuSig.LDA in directory ",outDir,"\n"))
        print("\n==================================\n\n")

        ## Run model fitting
        ## After the model is fitted, we may obtain extracted signatures and inferred exposures
        X = format_counts_lda(snv_counts);

        ## Save the likelihood for each K.
        likelihoods = NamedArray( repeat([NaN],9) , (string.([2:10;])) );
        for K in 2:10
            print("\nCalculating likelihood for K = "*string(K)*"...\n")
            ## Specify seed used before fitting model for each dataset,
            ## for each signature number (K)
            Random.seed!(seedInUse);
            model = LDA(K, 0.1, 0.1, 96, X);
            fit!(model, tol=1e-5);
            likelihood = MultiModalMuSig.calculate_loglikelihood(model);
            likelihoods[string(K)] = likelihood;
        end

        ## Elbow method: Find the inflection point of log-likelihood
        ## Calculate the first derivative of likelihoods
        ## Set names of deriv1 as "2"..."10"
        deriv1 = NamedArray( repeat([NaN],9) , (string.([2:10;])) );    
        for K in 2:10
          if K == 2
            ## For the smallest possible K specified by user,
            ## calculate 1st-derivative using forward difference operator
            ## with spacing equals to 1.
            deriv1[string(K)] = likelihoods[string(K+1)] - likelihoods[string(K)];
          elseif K == 10 
            ## For the largest possible K,
            ## calculate 1st-derivative using backward difference operator.
            deriv1[string(K)] = likelihoods[string(K)] - likelihoods[string(K-1)];
          else ## Calculate 1st-derivative using central difference
            deriv1[string(K)] = (likelihoods[string(K+1)] - likelihoods[string(K-1)]) / 2;
          end
        end

        ## Calculate the second derivative of likelihoods
        ## Set names of deriv2 as "2"..."10"
        deriv2 = NamedArray( repeat([NaN],9) , (string.([2:10;])) );
        for K in 2:10
          if K == 2
            ## For the smallest possible K specified by user,
            ## calculate 1st-derivative of the 1st-derivative using forward difference operator
            ## with spacing equals to 1.
            deriv2[string(K)] = deriv1[string(K+1)] - deriv1[string(K)];
          elseif K == 10 
            ## For the largest possible K,
            ## calculate 1st-derivative of the 1st-derivative using backward difference operator.
            deriv2[string(K)] = deriv1[string(K)] - deriv1[string(K-1)];
          else ## Calculate 1st-derivative using central difference
            deriv2[string(K)] = (deriv1[string(K+1)] - deriv1[string(K-1)]) / 2;
          end
        end

        ## Choose the smallest K where 2nd-derivative of likelihood function
        ## changes sign in its neighborhood.
        KBest = NaN; # Must declare outside of for loop
        
        for K in 2:10
          ## If deriv2["2"] and deriv2["3"] have opposite sign
          ## set KBest = 2
          if K == 2
            if sign(deriv2[string(K)]) * sign(deriv2[string(K+1)]) == -1
              KBest = K;
              break
            end
          ## If deriv2[string(K-1)] and deriv2[string(K+1)] have opposite sign
          ## set KBest = 2
          elseif K < 10
            if sign(deriv2[string(K-1)]) * sign(deriv2[string(K+1)]) == -1
              KBest = K;
              break
            end
          else
            KBest = K;
          end
        end


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
        outDir = string(topLevelFolder4Run,"/MultiModalMuSig.CTM.results/",datasetName,"/seed.",seedInUse)
        snv_counts = CSV.read(outDir*"/ground.truth.syn.catalog.tsv", delim='\t');

        print("\n\n==================================\n")
        print(string("Running MultiModalMuSig.CTM in directory ",outDir,"\n"))
        print("\n==================================\n\n")

        ## Run model fitting
        ## After the model is fitted, we may obtain extracted signatures and inferred exposures
        X = format_counts_mmctm(snv_counts);

        ## Save the likelihood for each K.
        likelihoods = NamedArray( repeat([NaN],9) , (string.([2:10;])) );
        for K in 2:10
            print("\nCalculating likelihood for K = "*string(K)*"...\n")
            ## Specify seed used before fitting model for each dataset,
            ## for each signature number (K)
            Random.seed!(seedInUse);
            model = MMCTM([K], [0.1], X);
            fit!(model, tol=1e-5);
            likelihood = MultiModalMuSig.calculate_loglikelihoods(model);
            likelihoods[string(K)] = likelihood[1];
        end
        print("\n\nFinished Calculating likelihood\n\n")


        ## Elbow method: Find the inflection point of log-likelihood
        ## Calculate the first derivative of likelihoods
        ## Set names of deriv1 as "2"..."10"
        deriv1 = NamedArray( repeat([NaN],9) , (string.([2:10;])) );    
        for K in 2:10
          if K == 2
            ## For the smallest possible K specified by user,
            ## calculate 1st-derivative using forward difference operator
            ## with spacing equals to 1.
            deriv1[string(K)] = likelihoods[string(K+1)] - likelihoods[string(K)];
          elseif K == 10 
            ## For the largest possible K,
            ## calculate 1st-derivative using backward difference operator.
            deriv1[string(K)] = likelihoods[string(K)] - likelihoods[string(K-1)];
          else ## Calculate 1st-derivative using central difference
            deriv1[string(K)] = (likelihoods[string(K+1)] - likelihoods[string(K-1)]) / 2;
          end
        end

        ## Calculate the second derivative of likelihoods
        ## Set names of deriv2 as "2"..."10"
        deriv2 = NamedArray( repeat([NaN],9) , (string.([2:10;])) );
        for K in 2:10
          if K == 2
            ## For the smallest possible K specified by user,
            ## calculate 1st-derivative of the 1st-derivative using forward difference operator
            ## with spacing equals to 1.
            deriv2[string(K)] = deriv1[string(K+1)] - deriv1[string(K)];
          elseif K == 10 
            ## For the largest possible K,
            ## calculate 1st-derivative of the 1st-derivative using backward difference operator.
            deriv2[string(K)] = deriv1[string(K)] - deriv1[string(K-1)];
          else ## Calculate 1st-derivative using central difference
            deriv2[string(K)] = (deriv1[string(K+1)] - deriv1[string(K-1)]) / 2;
          end
        end

        ## Choose the smallest K where 2nd-derivative of likelihood function
        ## changes sign in its neighborhood.
        KBest = NaN; # Must declare outside of for loop
        
        for K in 2:10
          ## If deriv2["2"] and deriv2["3"] have opposite sign
          ## set KBest = 2
          if K == 2
            if sign(deriv2[string(K)]) * sign(deriv2[string(K+1)]) == -1
              KBest = K;
              break
            end
          ## If deriv2[string(K-1)] and deriv2[string(K+1)] have opposite sign
          ## set KBest = 2
          elseif K < 10
            if sign(deriv2[string(K-1)]) * sign(deriv2[string(K+1)]) == -1
              KBest = K;
              break
            end
          else
            KBest = K;
          end
        end

        KBest = [KBest];


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
