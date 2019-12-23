#!/bin/bash
# get all filename in specified path

dataset=$1
l_path=$2
u_path=$3
rowperm_path=$4
colperm_path=$5

./sn_computing $dataset/ASIC_100k.mtx $l_path/ASIC_100k.l.sort.final $u_path/ASIC_100k.u.sort.final $rowperm_path/ASIC_100k.rowperm $colperm_path/ASIC_100k.colperm $6 $7
./sn_computing $dataset/ASIC_100ks.mtx $l_path/ASIC_100ks.l.sort.final $u_path/ASIC_100ks.u.sort.final $rowperm_path/ASIC_100ks.rowperm $colperm_path/ASIC_100ks.colperm $6 $7
./sn_computing $dataset/ASIC_320k.mtx $l_path/ASIC_320k.l.sort.final $u_path/ASIC_320k.u.sort.final $rowperm_path/ASIC_320k.rowperm $colperm_path/ASIC_320k.colperm $6 $7
./sn_computing $dataset/ASIC_320ks.mtx $l_path/ASIC_320ks.l.sort.final $u_path/ASIC_320ks.u.sort.final $rowperm_path/ASIC_320ks.rowperm $colperm_path/ASIC_320ks.colperm $6 $7
./sn_computing $dataset/ASIC_680k.mtx $l_path/ASIC_680k.l.sort.final $u_path/ASIC_680k.u.sort.final $rowperm_path/ASIC_680k.rowperm $colperm_path/ASIC_680k.colperm $6 $7
./sn_computing $dataset/ASIC_680ks.mtx $l_path/ASIC_680ks.l.sort.final $u_path/ASIC_680ks.u.sort.final $rowperm_path/ASIC_680ks.rowperm $colperm_path/ASIC_680ks.colperm $6 $7
./sn_computing $dataset/bcircuit.mtx $l_path/bcircuit.l.sort.final $u_path/bcircuit.u.sort.final $rowperm_path/bcircuit.rowperm $colperm_path/bcircuit.colperm $6 $7
./sn_computing $dataset/circuit_1.mtx $l_path/circuit1.l.sort.final $u_path/circuit1.u.sort.final $rowperm_path/circuit_1.rowperm $colperm_path/circuit_1.colperm $6 $7
./sn_computing $dataset/circuit_2.mtx $l_path/circuit2.l.sort.final $u_path/circuit2.u.sort.final $rowperm_path/circuit_2.rowperm $colperm_path/circuit_2.colperm $6 $7
./sn_computing $dataset/circuit_3.mtx $l_path/circuit3.l.sort.final $u_path/circuit3.u.sort.final $rowperm_path/circuit_3.rowperm $colperm_path/circuit_3.colperm $6 $7
./sn_computing $dataset/circuit_4.mtx $l_path/circuit4.l.sort.final $u_path/circuit4.u.sort.final $rowperm_path/circuit_4.rowperm $colperm_path/circuit_4.colperm $6 $7
./sn_computing $dataset/ckt11752_dc_1.mtx $l_path/cktdc.l.sort.final $u_path/cktdc.u.sort.final $rowperm_path/ckt11752_dc_1.rowperm $colperm_path/ckt11752_dc_1.colperm $6 $7
./sn_computing $dataset/ckt11752_tr_0.mtx $l_path/ckttr.l.sort.final $u_path/ckttr.u.sort.final $rowperm_path/ckt11752_tr_0.rowperm $colperm_path/ckt11752_tr_0.colperm $6 $7
./sn_computing $dataset/coupled.mtx $l_path/coupled.l.sort.final $u_path/coupled.u.sort.final $rowperm_path/coupled.rowperm $colperm_path/coupled.colperm $6 $7
./sn_computing $dataset/G2_circuit.mtx $l_path/g2.l.sort.final $u_path/g2.u.sort.final $rowperm_path/G2_circuit.rowperm $colperm_path/G2_circuit.colperm $6 $7
./sn_computing $dataset/G3_circuit.mtx $l_path/g3.l.sort.final $u_path/g3.u.sort.final $rowperm_path/G3_circuit.rowperm $colperm_path/G3_circuit.colperm $6 $7
./sn_computing $dataset/hcircuit.mtx $l_path/hcircuit.l.sort.final $u_path/hcircuit.u.sort.final $rowperm_path/hcircuit.rowperm $colperm_path/hcircuit.colperm $6 $7
./sn_computing $dataset/memplus.mtx $l_path/memplus.l.sort.final $u_path/memplus.u.sort.final $rowperm_path/memplus.rowperm $colperm_path/memplus.colperm $6 $7
./sn_computing $dataset/onetone1.mtx $l_path/onetone1.l.sort.final $u_path/onetone1.u.sort.final $rowperm_path/onetone1.rowperm $colperm_path/onetone1.colperm $6 $7
./sn_computing $dataset/onetone2.mtx $l_path/onetone2.l.sort.final $u_path/onetone2.u.sort.final $rowperm_path/onetone2.rowperm $colperm_path/onetone2.colperm $6 $7
./sn_computing $dataset/Raj1.mtx $l_path/raj1.l.sort.final $u_path/raj1.u.sort.final $rowperm_path/Raj1.rowperm $colperm_path/Raj1.colperm $6 $7
./sn_computing $dataset/rajat03.mtx $l_path/rajat03.l.sort.final $u_path/rajat03.u.sort.final $rowperm_path/rajat03.rowperm $colperm_path/rajat03.colperm $6 $7
./sn_computing $dataset/rajat15.mtx $l_path/rajat15.l.sort.final $u_path/rajat15.u.sort.final $rowperm_path/rajat15.rowperm $colperm_path/rajat15.colperm $6 $7
./sn_computing $dataset/rajat16.mtx $l_path/rajat16.l.sort.final $u_path/rajat16.u.sort.final $rowperm_path/rajat16.rowperm $colperm_path/rajat16.colperm $6 $7
./sn_computing $dataset/rajat17.mtx $l_path/rajat17.l.sort.final $u_path/rajat17.u.sort.final $rowperm_path/rajat17.rowperm $colperm_path/rajat17.colperm $6 $7
./sn_computing $dataset/rajat18.mtx $l_path/rajat18.l.sort.final $u_path/rajat18.u.sort.final $rowperm_path/rajat18.rowperm $colperm_path/rajat18.colperm $6 $7
./sn_computing $dataset/rajat20.mtx $l_path/rajat20.l.sort.final $u_path/rajat20.u.sort.final $rowperm_path/rajat20.rowperm $colperm_path/rajat20.colperm $6 $7
./sn_computing $dataset/rajat21.mtx $l_path/rajat21.l.sort.final $u_path/rajat21.u.sort.final $rowperm_path/rajat21.rowperm $colperm_path/rajat21.colperm $6 $7
./sn_computing $dataset/rajat22.mtx $l_path/rajat22.l.sort.final $u_path/rajat22.u.sort.final $rowperm_path/rajat22.rowperm $colperm_path/rajat22.colperm $6 $7
./sn_computing $dataset/rajat23.mtx $l_path/rajat23.l.sort.final $u_path/rajat23.u.sort.final $rowperm_path/rajat23.rowperm $colperm_path/rajat23.colperm $6 $7
./sn_computing $dataset/rajat24.mtx $l_path/rajat24.l.sort.final $u_path/rajat24.u.sort.final $rowperm_path/rajat24.rowperm $colperm_path/rajat24.colperm $6 $7
./sn_computing $dataset/rajat25.mtx $l_path/rajat25.l.sort.final $u_path/rajat25.u.sort.final $rowperm_path/rajat25.rowperm $colperm_path/rajat25.colperm $6 $7
./sn_computing $dataset/rajat26.mtx $l_path/rajat26.l.sort.final $u_path/rajat26.u.sort.final $rowperm_path/rajat26.rowperm $colperm_path/rajat26.colperm $6 $7
./sn_computing $dataset/rajat27.mtx $l_path/rajat27.l.sort.final $u_path/rajat27.u.sort.final $rowperm_path/rajat27.rowperm $colperm_path/rajat27.colperm $6 $7
./sn_computing $dataset/rajat28.mtx $l_path/rajat28.l.sort.final $u_path/rajat28.u.sort.final $rowperm_path/rajat28.rowperm $colperm_path/rajat28.colperm $6 $7
./sn_computing $dataset/rajat29.mtx $l_path/rajat29.l.sort.final $u_path/rajat29.u.sort.final $rowperm_path/rajat29.rowperm $colperm_path/rajat29.colperm $6 $7

