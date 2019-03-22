#!/bin/sh

set -e

dir="mc"

nev=2000000
numu_xse="$GENIE/share/xsec/numu_C12.xml"
anmu_xse="$GENIE/share/xsec/anmu_C12.xml"
seed=25
messenger=""
# messenger="$PWD/Messenger.xml"
redir="/dev/null"

ty="numu"
# numu-C12 1 GeV
mkdir -p $dir/ev0_${ty}_1GeV && cd $dir/ev0_${ty}_1GeV
echo "gevgen -r 0 -n $nev -e 1 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 0 -n $nev -e 1 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# numu-C12 3 GeV
mkdir -p $dir/ev1_${ty}_3GeV && cd $dir/ev1_${ty}_3GeV
echo "gevgen -r 1 -n $nev -e 3 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 1 -n $nev -e 3 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# numu-C12 50 GeV
mkdir -p $dir/ev2_${ty}_50GeV && cd $dir/ev2_${ty}_50GeV
echo "gevgen -r 2 -n $nev -e 50 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 2 -n $nev -e 50 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# numu-C12 80 GeV
mkdir -p $dir/ev3_${ty}_80GeV && cd $dir/ev3_${ty}_80GeV
echo "gevgen -r 3 -n $nev -e 80 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 3 -n $nev -e 80 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# numu-C12 100 GeV
mkdir -p $dir/ev4_${ty}_100GeV && cd $dir/ev4_${ty}_100GeV
echo "gevgen -r 4 -n $nev -e 100 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 4 -n $nev -e 100 -p 14 -t 1000060120 --cross-sections $numu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

ty="anmu"
# anmu-C12 1 GeV
mkdir -p $dir/ev0_${ty}_1GeV && cd $dir/ev0_${ty}_1GeV
echo "gevgen -r 0 -n $nev -e 1 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 0 -n $nev -e 1 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# anmu-C12 3 GeV
mkdir -p $dir/ev1_${ty}_3GeV && cd $dir/ev1_${ty}_3GeV
echo "gevgen -r 1 -n $nev -e 3 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 1 -n $nev -e 3 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# anmu-C12 50 GeV
mkdir -p $dir/ev2_${ty}_50GeV && cd $dir/ev2_${ty}_50GeV
echo "gevgen -r 2 -n $nev -e 50 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 2 -n $nev -e 50 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# anmu-C12 80 GeV
mkdir -p $dir/ev3_${ty}_80GeV && cd $dir/ev3_${ty}_80GeV
echo "gevgen -r 3 -n $nev -e 80 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 3 -n $nev -e 80 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -

# anmu-C12 100 GeV
mkdir -p $dir/ev4_${ty}_100GeV && cd $dir/ev4_${ty}_100GeV
echo "gevgen -r 4 -n $nev -e 100 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger"
nohup gevgen -r 4 -n $nev -e 100 -p -14 -t 1000060120 --cross-sections $anmu_xse --seed $seed --message-thresholds $messenger > $redir &
cd -
