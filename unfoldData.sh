#/usr/bin/sh
python Unfolding.py \
    legotrain_610_20181010-1926_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_610_20181010-1926_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_CF_pPb_2305_20190109_LHC13bcde.root \
    MB

python Unfolding.py \
    legotrain_610_20181010-1926_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_610_20181010-1926_LHCb4_fix_CF_pPb_MC_ptHardMerged.root \
    legotrain_CF_pPb_2305_20190109_LHC13bcde.root \
    Triggered
