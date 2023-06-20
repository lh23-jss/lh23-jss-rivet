#!/bin/bash

set -x

echo "args: $@"

suffix=

# PN, ParT
model=$1
if [[ "$model" == "ParT" ]]; then
    modelopts="example_ParticleTransformer.py --use-amp --optimizer-option weight_decay 0.01"
    lr="1e-3"
elif [[ "$model" == "PN" ]]; then
    modelopts="example_ParticleNet.py"
    lr="1e-2"
else
    echo "Invalid model $model!"
    exit 1
fi

weaver \
    --data-train "W:${DATADIR}/ZW_nunuqq/WZ_*.root" "QCD:${DATADIR}/Dijet_Slice_600upGeV/DIJET_*.root" \
    --data-test "W:${DATADIR}/ZW_nunuqq/WZ_*.root" "QCD:${DATADIR}/Dijet_Slice_600upGeV/DIJET_*.root" \
    --data-config WvsQCD.yaml --network-config $modelopts \
    --model-prefix training/WvsQCD/${model}/{auto}${suffix}/net \
    --num-workers 1 --fetch-step 1 --in-memory --train-val-split 0.8 \
    --batch-size 512 --samples-per-epoch 512000 --samples-per-epoch-val 128000 --num-epochs 20 --gpus 0 \
    --start-lr $lr --optimizer ranger --log logs/WvsQCD_${model}_{auto}${suffix}.log --predict-output pred.root \
    --cross-validation 'idx%5' \
    "${@:2}"
