#/usr/bin/sh
python Unfolding.py \
    legotrain_512_20180523-2331_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_512_20180523-2331_LHCb4_fix_CF_pPb_MC_ptHardMerged.root\
    legotrain_CF_pPb_1839_20180613_LHC13bcde.root \
    MB

python Unfolding.py \
    legotrain_512_20180523-2331_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_512_20180523-2331_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_CF_pPb_1839_20180613_LHC13bcde.root \
    Triggered
