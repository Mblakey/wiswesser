/*********************************************************************
 
Author : Michael Blakey

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "parser.h"

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

const char *cli_inp;
const char *format; 

#define BENCH_N 417

/* not fast, but will do the job */ 
void RunBenchmark() {
  
  const char *wln_bench[BENCH_N]  = 
  {
    "1V1",
    "2O2",
    "Q2Q",
    "WN1O1NW",
    "Z3Z",
    "QV3VQ",
    "ZV2VZ",
    "ZVM1MVZ",
    "3O2",
    "QV1",
    "2M1",
    "Z2VQ",
    "2VO2",
    "1SVM1",
    "1VOV1",
    "Q2G",
    "3M2",
    "Z2E",
    "QV4",
    "E2E",
    "ZQ",
    "ZH",
    "QH",
    "IH",
    "QV1G",
    "ZV2M1VQ",
    "2OVM1VO2",
    "12VM1",
    "G1OS2S2",
    "QNU3",
    "2U2U1",
    "MU2O1",
    "ZVMNU2",
    "SU3",
    "1UU1",
    "3UU2",
    "G1UU2",
    "G1UU1G",
    "1NU1UN1",
    "1U1U1",
    "OCO",
    "NCH",
    "QNO",
    "WNQ",
    "ONONO",
    "ZNU1",
    "2CN",
    "2OCN",
    "2NCO",
    "SCN3",
    "NCS3",
    "OC2",
    "GH",
    "EH",
    "2H",
    "10H",
    "VHVH",
    "VH3",
    "SH2Q",
    "SH2G",
    "SH1M2",
    "2Y2",
    "GXGGG",
    "G2N2&3",
    "2OPO&2&O2",
    "Q1XGG2Y1Q1Z",
    "1S2YZVMYVQS1",
    "W-G-QO",
    "F-G-FF",
    "12SW12",
    "WS2&12",
    "QVY9&19",
    "12M1",
    "12N3&2",
    "19YQM1",
    "QVY19&2Q",
    "G1XGGYP3&3&&P2&2",
    "OK2&2&2",
    "WSQO1",
    "2B2&2",
    "ZXZUS",
    "ZSW2",
    "QVYZY2&2",
    "QY",
    "1Y&N1&1",
    "2OVX",
    "QY&1X&&1Y",
    "Z2X2&&2Y2",
    "GXGG2NY&&Y",
    "2Y2&3YX&&&1VOY",
    "1Y&N1VM1O1O1&Y",
    "OV1 &-NA-",
    "SUXS&O4 &-KA-",
    "ZVSH &ZV1Z",
    "Z6Z &Q2 &Q2",
    "Z3Z &QH &QH",
    "QV2 &ZH",
    "WSQ2G &ZH",
    "3M3 &WSQQ",
    "4N4&4 &GH",
    "Z3Z &GH &GH &QH &QH",
    "QVVQ &ZH &ZH",
    "1K &G",
    "G2SH2G2G &G",
    "G2KO&2&2G &GH",
    "OV1 &OV1 &-CA-",
    "1R",
    "GV1N2&R",
    "QY2&XQR&R",
    "RM1R",
    "WNR CVM1VQ",
    "Z1YQR DO1",
    "QVR BNUNR DN1&1",
    "2N2&R COVMR DO1",
    "WNR DMNU1R CNW",
    "WNR BMNU1R CNW",
    "WNR CMNU1R BNW",
    "WNR DMNU1R DNW",
    "ZSWR DSZW",
    "ZMR CG",
    "QMR BE",
    "VHR CSH",
    "SHR DVQ",
    "QV1N1VQR BVQ",
    "ZSWR CN1O2&1O2",
    "2O1N1O2&R B2O1",
    "ZSWR D-AS-U-AS-R DSWQ",
    "ZR BG DE",
    "ZR CG BE",
    "ZR BQ ENW",
    "ZR BG EVQ",
    "WNR CG FE DOV1",
    "3OR BF E2 CM1",
    "1VOR CG EF BO1 DVO1",
    "WNR BQ CNW ENW",
    "QVR BG FVQ",
    "WNR CNW DMNU2",
    "ZR C1U1R DZ",
    "WNR CF DI E2U1",
    "QMR DQ CO1 EO1",
    "Q-AS-QO&R BF CF EF",
    "ZR BG EG CE FE",
    "SHR DE CV1",
    "VHR DG CVG E1G",
    "WSQR BO2 ESWQ",
    "1N1&R CYR&R DN1&1",
    "GR DG BOR BO1 EO1",
    "WNR DNW BR BG CNW",
    "G1O1VR CO1O1G E1R CG DG",
    "GXGGR B1O1O1Q DXGGG",
    "WNR BR& ENW",
    "1VOR BVR& EOV1",
    "G1OVR BR CQ& DVO1",
    "WNR DNW BR CQ",
    "RR BR BR DR&& ER CR",
    "ZR BYR DZ&1NV1&V1",
    "WNR BG EYR DZ CG&1N2Q2Q",
    "G1Y&M1VR CV1MR DG& DN1&1",
    "T6NJ",
    "L6V DVJ",
    "L6VVJ",
    "T5M CNJ",
    "T6ON DMJ",
    "T6NSO ENJ",
    "T6NSN EMJ",
    "T5NO DNJ",
    "T7M CN ENJ",
    "T5NS DS CHJ",
    "T6M CM DMTJ",
    "T5N CO EUTJ",
    "L6U CUTJ",
    "L5 AHJ",
    "T5VOVTJ",
    "T5NN CHJ",
    "T6OV EUTJ",
    "T6NJ CVQ DG",
    "T5OJ BVH DVQ",
    "L5TJ AVQ AVQ B1 B1 CVQ",
    "T6NJ BVQ FG",
    "T5MVMVTJ EVM1 EQ",
    "T6VMVMVJ F2 F2",
    "T5UTJ A1 B1 C1 C1",
    "L5UTJ A2 DVQ E1 E1",
    "T5NTJ A1 C1 DR",
    "L6Y DYJ AUNQ DUNQ",
    "L6Y BUTJ AUNMR& EG FG",
    "T5SWTJ",
    "T5-AS-J",
    "T6 DUJ",
    "T6P DUJ",
    "T6P-KA- DUJ",
    "T6S-SB-HS DHJ",
    "T-10-M FOTJ",
    "L6V BUTJ",
    "L6V BUTJ EG FG",
    "L5UTJ A1 DV1Q E1 E1",
    "L5TJ AVQ A1 B1 B1 CVQ",
    "T5NYMVJ A1 BUM E1R",
    "T5NYVOJ BU1Q E1R",
    "T5NNVYJ DUNR& E1",
    "T6VMVNVJ D1 FMVZ",
    "T6VOV EOTJ DR& FR",
    "T6OPOTJ",
    "T6-SI- C-SI- E-SI-TJ",
    "L3TJ AO1",
    "T5NNV AHJ A1 BR DE& E1",
    "T5M CN BUTJ B1NR&1R",
    "L7UUTJ",
    "L6V BUTJ B1 EYU1",
    "L6UTJ C1 C1 DVQ D1",
    "L6V DYJ BVQ DUNZ",
    "T5OSOTJ BO",
    "T5NNV EHJ BNO D1 E1",
    "T5NYMV EHJ A1 BUM",
    "T6OYOYOYTJ BUM DUM FUM",
    "T6N CNJ BZ DZ ER DG& F2",
    "T6VMVMV FHJ F2 FR",
    "T6NJ C2Q E1Q &GH",
    "T6KJ A1R& COVN1&1 &E",
    "T6N DNTJ AYR&R DG& D1 &GH &GH",
    "T5K CSJ A1R BZ& E1 &GH &G",
    "T6KJ AO",
    "L66TJ",
    "L66J",
    "T56 CMJ",
    "T56 BOJ",
    "T55 AN DM FST&J &G",
    "L56 BHJ",
    "L66J CVQ",
    "T56 BMJ HO1",
    "L66J CVQ DG",
    "L56 CHJ B1 C1 C1",
    "T66 BN GNJ EVQ",
    "T66 BN&TJ",
    "T66 BN BU DUTJ",
    "T56 BN CH&TJ",
    "T66 CM AU- FTJ",
    "T66 BN ES DHJ",
    "T66 BS DN BHJ",
    "T55 AK DM FST&J &G",
    "T55 AN CM GNT&J",
    "T66 BSN EM&TJ",
    "T6-10- IM&TJ",
    "T55 CN GNJ BE",
    "L66&TJ",
    "L66J BQ EQ",
    "L66&TJ GQ JG",
    "T56 BO AU- E GUTJ",
    "L B666J",
    "L C666J",
    "L D6 C6 B566 MNJ",
    "T B566 DMJ",
    "T E6 C666 ANV&&TTJ",
    "T C566 DMJ",
    "L B656 HVJ ENW KG",
    "T I5 F6 D6 C665 JO VOJ",
    "T C666 BK ISJ B2N1&1",
    "T C666 BNJ B2N1&1",
    "T C666 BKJ B1 EZ MZ &GH &G",
    "L B5 G566J",
    "T D6 B666 KN&TT&J",
    "L E5 B666 OV MUTJ A1 E1",
    "L C666 BV IVJ D1Q GQ KQ",
    "L F6 E6 B666 BUTJ F1 IQ J1 J1 N1 O1 R1 U1 U1",
    "T C6 B566 LO&TT&J EQ FQ JQ OQ",
    "L566 1A LJ",
    "L666 B6 2AB PJ",
    "L C6665 1A PJ",
    "T666 B6 C6 3ABC S KMJ",
    "L C6566 O6 U6 N5 2AO C&J",
    "L C6666 P6 W6 O6 2AP E&-&&&TTTTJ",
    "L666 N6 B6 R6 Q6 4ABRS C&J",
    "L F6 D6 B6 P6666 2AB E&J",
    "T G6 E6 B6666 C6 3ABC D& FN MN WHJ",
    "L D66 K666 1A U EHJ",
    "T D6 C6666 B6 T6 5ABCDU B& MOJ",
    "L E6 C6666 T6 S6 2AT E&J",
    "T F6 C6 B5 O6666 D6 4ABCD F& MSJ",
    "L E6 E-6 D5 C4566 2AE- WTJ",
    "L666 B-6 B6 C6 4ABB-C U GH MHJ",
    "L666 B-6 B6 B-&6 B6 C6 6ABBB-B-&C D& BXTJ",
    "L665 K6 J6 J-6 I6 H6 6AIJJJ-K D& JXTJ",
    "T E6 E-6 C6666 2AE- W CN ON CHJ",
    "T E6 E-6 D5 C5566 2AE- A& OSJ",
    "T656 M6 B6 Q6 P5 4ABQR A& HN WNJ",
    "L H6 H--6 H-6 H--&6 G6 E6 C6666 B-6 B6 7ABB-CH-H--H--& P&-TJ",
    "T66 BNJ EQ HO1 IQ D- CT6NJ EVQ",
    "T6N DNJ B- CT6NTJ A1 EQ",
    "T6NTJ A1 B- BT6NJ",
    "L66J C- DL66J B- AL6TJ",
    "T66 CNJ H- HT66 CNJ G- AL6TJ C- AL6TJ",
    "T5SJ B- CT5SJ",
    "T6NJ BN2N1&1&1- BT5SJ EE",
    "T6NJ BO2N2&2 FXQR DR&&- BT5OJ",
    "L66J C- DL6OJ& D- DT6NJ",
    "L66J C- DT6OTJ& D- DT6NJ",
    "L66J BG EE C- DT6NJ& D- DT6OTJ",
    "T6-GE- DOTJ A-& AT6-GE-TJ",
    "T66 COJ EVQ H-& DT6OJ",
    "T6SXS EXTJ B-& AL5XTJ& E-& AL5XTJ",
    "T56 BO-AS-OJ C-& CT56 BO-AS-OJ",
    "T E6 D6 B666 CO NXJ N-& BT56 BXO DHJ",
    "T5VOXJ C-& AL4XXTJ",
    "T666 1A N AX DMTJ",
    "L55 ATJ",
    "T55 COTJ",
    "T55 A COTJ",
    "L66 A BTJ",
    "T66 A B AMTJ",
    "T56 A ANTJ A1 GQ",
    "T67 A B AO HOTJ",
    "T66 A DM HOTJ",
    "T66 A B ASTJ",
    "T56 A DOTJ",
    "T56 A FOTJ",
    "T67 A DMTJ",
    "L666/GL 2AF LTJ",
    "L556/EK 2AE K EXTJ",
    "L5 G66 N66/FS 2AE STJ",
    "L646/B-F/BI A A 2BF I BX FXTJ",
    "L D55 D6 2DH L DX HXTJ",
    "L566 F6/FK/FO 3AEF O FXTJ",
    "L535 B5/CG 3ABC ITJ",
    "T C655 A ASJ",
    "T C655 A AM EO GO FHJ",
    "L D6 C555 A AH JHJ",
    "T C555 B6 A 1B N AOX EOJ",
    "T57 B5 A 1B L BXOTJ",
    "T C6566 N6 O 1A S IN QN AU PUTJ",
    "T D3556 J5 K 2AF N EOX HO NO GH KHJ",
    "T666 B 1A L DOTJ",
    "T76 B6 A C 1B L EO FHTJ",
    "L666/GL 2AF LTJ",
    "L C66 L666/KT 2AJ TTJ",
    "T C66 L56 R66/KW 2AJ W NOTJ",
    "T D666 C6 A B 1C P CX JU NH&TT&J",
    "T6 G5 F66 C6 A B 2CF S AOOX FX IOTJ",
    "T566 F5/FK/FN 3AEF N CO FXTJ",
    "T5536/FK 3AAE K AX COTJ",
    "L6 H665/FP 2AF P FXTJ",
    "T54 H5 B6 B6/GO A 3BBF O BX GN AH&&&&TJ",
    "T B666/GL A 2BG L GXOO JHTJ",
    "L F5 E5 C655 A F-TJ",
    "L E5 D5 C555 A E-TJ",
    "L F5 E5 D5 C555 A D- F-TJ",
    "L55 F6 F6/BG/HM A 2FG MTJ",
    "L C656 M6 J5 A KTJ",
    "L E6 D58 N65 A D-TJ",
    "L C657 N6 K5 A LTJ",
    "L66 B6 A B- C 1B ITJ",
    "T66 B6 A B- C 1B I BN DN FN HNTJ",
    "T66 B6 A B- C 1B ITJ BE",
    "L555/BG 2AB ITJ",
    "T566/BH 2AB K AN JMTJ",
    "L545/BF 2AB HTJ",
    "T555/BG 2AB I DOTJ",
    "T C3666/FM 2AG M DO LOTJ",
    "T B5565/GM 2AH M IO KOTJ",
    "T666 B C 1A K EOTJ",
    "L666 B C 1A K EHTJ",
    "T756/EK F 2AG K HOTJ",
    "T B4656/FL 3ABE L IOTJ",
    "T B66 L667/HT 3ABG T KM RKHOSTJ",
    "T C6 B666/JO A 2BK O AOTJ",
    "T76 B6 A C 1B L EOTJ",
    "T D5 B665 A 1B N EO GO MNOTJ",
    "T E6 D566 B5 A C 2BD Q AM DX HMTJ",
    "T B6 B6566 E5/FR 5ABBFG T BX EN LOTJ",
    "T-11-56/IO J 2AK O ANOSO HO MHTJ",
    "L435 B3 2AB GTJ",
    "T655 C6 B D 2AC L BOTJ",
    "T C55 F555 B5 D I 5ABCGH M LOTJ",
    "L B4 F6 D7 C4 L676 M P 4CDNO VTJ",
    "L55 F56 B5 C G I 3ABH LTJ",
    "T665 C5/EK D 3AEF L DO GO JNTJ",
    "T C5 B5535/HL 3ABI L EOTJ",
    "L455 B5 C5 E5 F 6ABCDEG LTJ",
    "L B5 E3555 C 3ABF KTJ",
    "L545 B4 C5 D 4ABCE JTJ",
    "L454 B5 C5 E 4ABCD JTJ",
    "T B555 C5 A D 2BC J HNTJ",
    "T C555 B5 D5 A E 3BCD LTJ",
    "T C5 B5535/HL 3ABI L CN EM GNTJ",
    "T555 B6/CI 3ABC L FO KOTJ",
    "T B666 B6/CL/HM A 4BBCG M AOTJ",
    "L B5555 D5/GL E 4AFGH LTJ",
    "T656 C5/DJ B 3ACD L GO LOTJ",
    "L E3 D5 D5 C555/FJ/BN A 3DEI NTJ",
    "T5 F6 E56 B6 B6/CR/NS A 5BBCEO S EX LM ONTJ",
    "T-T C666 DM ISJ IH IOSR D1UN- F-12-J",
    "T-T E6 C66 P56 V66/OD& 2AN D& ROTJ K- CT5OTJ D- D-6-TJ",
    "T-L C666 DY GYJ GUN- DL C666J GNU- D-10-J",
    "T-T56 CMJ D1- BT56 CMJ D1- BT56 CMJ D1- BT56 CMJ D1- B-16-J",
    "L G6 C-10-66&T&&J",
    "L H6 F6-11-6 ATJ",
    "T6 I6 M6-12-6/BU 4ABMN A& GSS SSSJ",
    "T-T56 BYN DHJ CR DR DNU- B-11-TJ",
    "T-L G6 E6 C666&&T&&J OOR DO- D-10-J",
    "T-T66 CNJ IO- JT66 CNJ B1R COR D1- B-18-J",
    "T-T55 BM HUTJ H- ET5NY EUTJ BU1- BT5MJ E1U- AT5YMYJ CU1- C-16-J",
    "L56 F0J 0-FE-- 0L56 B0TJ",
    "L56 F0J 0-FE-- 0L5 B0J A1",
    "T-L6 B0J A- AL5 B0J 0-FE-- 0L5 B0J A- AL6 B0J 0-FE-- 0-6-J",
    "T-L5 B0J A- AL5 B0J CG 0-FE-- 0L5 B0J A- AL5 B0J C1 0-FE-- 0-6-J",
    "D5ZD-CU-DZTJ &Q &Q",
    "D6O-AU-DVJ B1 B1 D1 F1",
    "D5ZD-CU-DZTJ &Q &Q",
    "D566 1A L BND-ZN-OJ C-& CD566 1A L BND-ZN-OJ",
    "D C65 J656 1A S AN HND-PT-DN IJ IG &G",
    "D B656 GND-FE-DNJ &WSW",
    "OV1 &-NA- &7/1",
    "T6KTJ A1 A1 &3/11",
    "T6KTJ A1 A1 &G &3/14",
    "T56 BS DSJ &G &6/13",
    "T56 CKT&J C-& AT6KTJ &E &6/23",
    "OVVO &-ZN- &8/1 &8/4",
    "T65 BNKJ C-& CT6TJ"
  }; 

  const char *smiles_bench[BENCH_N] = 
  {
    "CC(=O)C",
    "CCOCC",
    "OCCO",
    "[O-][N+](=O)COC[N+](=O)[O-]",
    "NCCCN",
    "OC(=O)CCCC(=O)O",
    "NC(=O)CCC(=O)N",
    "NC(=O)NCNC(=O)N",
    "CCOCCC",
    "CC(=O)O",
    "CNCC",
    "NCCC(=O)O",
    "CCOC(=O)CC",
    "CNC(=O)SC",
    "CC(=O)OC(=O)C",
    "OCCCl",
    "CCNCCC",
    "NCCBr",
    "CCCCC(=O)O",
    "BrCCBr",
    "NO",
    "N",
    "O",
    "I",
    "OC(=O)CCl",
    "NC(=O)CCNCC(=O)O",
    "CCOC(=O)NCC(=O)OCC",
    "CCCCCCCCCCCCC(=O)NC",
    "ClCOSCCSCC",
    "CCC=NO",
    "CC=CC=C",
    "COCC=N",
    "CC=NNC(=O)N",
    "CCC=S",
    "C#C",
    "CCC#CC",
    "CC#CCl",
    "ClC#CCl",
    "CN=C=NC",
    "C=C=C",
    "O=C=O",
    "C#N",
    "ON=O",
    "[O-][N+](=O)O",
    "O=NON=O",
    "NN=C",
    "CCC#N",
    "CCOC#N",
    "CCN=C=O",
    "CCCN=C=S",
    "CCCSC#N",
    "CC=C=O",
    "Cl",
    "Br",
    "CC",
    "CCCCCCCCCC",
    "O=CC=O",
    "CCCC=O",
    "OCCS",
    "SCCCl",
    "CCNCS",
    "CCC(CC)C",
    "ClC(Cl)(Cl)Cl",
    "CCCN(CCCl)CC",
    "CCP(=O)(OCC)OCC",
    "NCC(CCC(CO)(Cl)Cl)CO",
    "CSCCC(C(=O)NC(C(=O)O)SC)N",
    "[O-][Cl](=O)(=O)O",
    "F[Cl](F)F",
    "CCCCCCCCCCCCS(=O)(=O)CCCCCCCCCCCC",
    "CCCCCCCCCCCCS(=O)(=O)CC",
    "CCCCCCCCCCCCCCCCCCCC(C(=O)O)CCCCCCCCC",
    "CCCCCCCCCCCCNC",
    "CCCCCCCCCCCCN(CCC)CC",
    "CCCCCCCCCCCCCCCCCCCC(NC)O",
    "CCCCCCCCCCCCCCCCCCCC(C(=O)O)CCO",
    "CCCP(C(C(CCl)(Cl)Cl)P(CC)CC)CCC",
    "CC[N+](CC)(CC)[O-]",
    "COS(=O)(=O)O",
    "CCB(CC)CC",
    "NC(=S)N",
    "CCS(=O)(=O)N",
    "CCC(C(C(=O)O)N)CC",
    "CC(O)C",
    "CN(C(C)C)C",
    "CCOC(=O)C(C)(C)C",
    "CC(CC(CC(O)C)(C)C)C",
    "NCCC(CCC(CC)C)(CC)C",
    "CC(N(C(C)C)CCC(Cl)(Cl)Cl)C",
    "CCC(CCCC(C(C)(C)C)CC(=O)OC(C)C)CC",
    "COCOCNC(=O)CN(C(C)C)C(C)C",
    "[O-]C(=O)C.[Na]",
    "CCCCOC(=S)[S].[K]",
    "NC(=O)S.NCC(=O)N",
    "NCCCCCCN.CCO.CCO",
    "NCCCN.O.O",
    "CCC(=O)O.N",
    "ClCCS(=O)(=O)O.N",
    "OS(=O)(=O)O.CCCNCCC",
    "CCCCN(CCCC)CCCC.Cl",
    "NCCCN.Cl.Cl.O.O",
    "OC(=O)C(=O)O.N.N",
    "C[N+](C)(C)C.[Cl-]",
    "ClCCS(CCCl)CCCl.[Cl-]",
    "ClCC[N+](CCCl)(CC)[O-].Cl",
    "[O-]C(=O)C.[O-]C(=O)C.[Ca]",
    "Cc1ccccc1",
    "CCN(c1ccccc1)CC(=O)Cl",
    "CCC(C(c1ccccc1)(c1ccccc1)O)O",
    "c1ccc(cc1)NCc1ccccc1",
    "OC(=O)CNC(=O)c1cccc(c1)[N+](=O)[O-]",
    "NCC(c1ccc(cc1)OC)O",
    "CN(c1ccc(cc1)N=Nc1ccccc1C(=O)O)C",
    "COc1ccc(cc1)NC(=O)Oc1cccc(c1)N(CC)CC",
    "[O-][N+](=O)c1ccc(cc1)NN=Cc1cccc(c1)[N+](=O)[O-"	,
    "[O-][N+](=O)c1cccc(c1)C=NNc1ccccc1[N+](=O)[O-]",
    "[O-][N+](=O)c1cccc(c1)NN=Cc1ccccc1[N+](=O)[O-]",
    "[O-][N+](=O)c1ccc(cc1)NN=Cc1ccc(cc1)[N+](=O)[O-]",
    "NS(=O)(=O)c1ccc(cc1)S(=O)(=O)N",
    "NNc1cccc(c1)Cl",
    "ONc1ccccc1Br",
    "O=Cc1cccc(c1)S",
    "OC(=O)c1ccc(cc1)S",
    "OC(=O)CN(c1ccccc1C(=O)O)CC(=O)O",
    "CCOCN(c1cccc(c1)S(=O)(=O)N)COCC",
    "COCCc1ccccc1N(COCC)COCC",
    "NS(=O)(=O)c1ccc(cc1)[As]=[As]c1ccc(cc1)S(=O)(=O)O",
    "Brc1ccc(c(c1)Cl)N",
    "Brc1c(N)cccc1Cl",
    "[O-][N+](=O)c1ccc(c(c1)N)O",
    "OC(=O)c1ccc(c(c1)N)Cl",
    "CC(=O)Oc1cc(Br)c(cc1Cl)[N+](=O)[O-]",
    "CCCOc1cc(CC)cc(c1F)NC",
    "COC(=O)c1c(F)cc(c(c1Cl)OC)OC(=O)C",
    "[O-][N+](=O)c1cc([N+](=O)[O-])c(c(c1)[N+](=O)[O-])O",
    "OC(=O)c1c(Cl)cccc1C(=O)O",
    "CC=NNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]",
    "Nc1ccc(cc1)C=Cc1cccc(c1)N",
    "C=CCc1cc(cc(c1I)F)[N+](=O)[O-]",
    "COc1cc(NO)cc(c1O)OC",
    "Fc1cc(F)c(c(c1)[As](=O)(O)O)F",
    "Brc1cc(Cl)c(c(c1Cl)N)Br",
    "Sc1ccc(c(c1)C(=O)C)Br",
    "O=Cc1cc(CCl)c(c(c1)C(=O)Cl)Cl",
    "CCOc1ccc(cc1S(=O)(=O)O)S(=O)(=O)O",
    "CN(c1ccc(cc1)C(c1cccc(c1)N(C)C)c1ccccc1)C",
    "COc1ccc(cc1Oc1cc(Cl)ccc1Cl)OC",
    "[O-][N+](=O)c1ccc(cc1c1cccc(c1Cl)[N+](=O)[O-])[N+](=O)[O-]",
    "ClCOCOc1cc(Cc2ccc(c(c2)Cl)Cl)cc(c1)C(=O)COCCl",
    "OCOCOCc1cc(ccc1C(Cl)(Cl)Cl)C(Cl)(Cl)Cl",
    "[O-][N+](=O)c1cc(ccc1c1ccccc1)[N+](=O)[O-]",
    "CC(=O)Oc1cc(ccc1C(=O)c1ccccc1)OC(=O)C",
    "ClCOC(=O)c1ccc(cc1c1cccc(c1)O)C(=O)OC",
    "Oc1cccc(c1)c1cc(ccc1[N+](=O)[O-])[N+](=O)[O-]",
    "c1ccc(cc1)c1ccccc1c1cc(ccc1c1ccc(cc1)c1ccccc1)c1cccc(c1)c1ccccc1",
    "Nc1ccc(cc1)C(c1ccccc1N)CN(C(=O)C)C(=O)C",
    "OCCN(CC(c1ccc(c(c1)[N+](=O)[O-])Cl)c1ccc(c(c1)Cl)N)CCO",
    "ClCC(NCC(=O)c1ccc(c(c1)C(=O)CNc1ccc(cc1)Cl)N(C)C)C",
    "c1cccnc1",
    "O=C1C=CC(=O)C=C1",
    "O=C1C=CC=CC1=O",
    "c1ncc[nH]1",
    "O1C=CNC=N1",
    "S1OC=NC=N1",
    "N1=CNC=NS1",
    "c1ncno1",
    "N1=CN=CNC=C1",
    "S1CSC=N1",
    "N1CCNNC1",
    "C1OCC=N1",
    "C1CC=CC=C1",
    "C1C=CC=C1",
    "O=C1CCC(=O)O1",
    "C1C=CN=N1",
    "O=C1CCC=CO1",
    "OC(=O)c1cnccc1Cl",
    "O=Cc1occ(c1)C(=O)O",
    "OC(=O)C1CCC(C1(C)C)(C(=O)O)C(=O)O",
    "Clc1cccc(n1)C(=O)O",
    "CNC(=O)C1(O)NC(=O)NC1=O",
    "CCC1(CC)C(=O)NC(=O)NC1=O",
    "CC1=C(C)C(CC1)(C)C",
    "CCC1=CCC(C1(C)C)C(=O)O",
    "CC1CN(CC1c1ccccc1)C",
    "ON=C1C=CC(=NO)C=C1",
    "ClC1C(Cl)CC=CC1=NNc1ccccc1",
    "O=S1(=O)CCCC1",
    "C1C=CC=[As]1",
    "C1=CC=C=C=C1",
    "P1=CC=C=C=C1",
    "[K]1=PC=C=C=C1",
    "C1S[Sb]SC=C1",
    "N1CCCCOCCCC1",
    "O=C1CCCC=C1",
    "ClC1C(Cl)CC=CC1=O",
    "OCC(=O)C1CC=C(C1(C)C)C",
    "OC(=O)C1CCC(C1(C)C)(C)C(=O)O",
    "CN1C(=N)NC(=O)C1Cc1ccccc1",
    "OC=C1N=C(OC1=O)Cc1ccccc1",
    "CC1N=NC(=O)C1=Nc1ccccc1",
    "NC(=O)NC1C(=O)NC(=O)N(C1=O)C",
    "O=C1OC(=O)C(OC1c1ccccc1)c1ccccc1",
    "C1CO[P]OC1",
    "[Si]1C[Si]C[Si]C1",
    "COC1CC1",
    "Brc1ccc(cc1)n1c(=O)cc(n1C)C",
    "c1ccc(cc1)CN(c1ccccc1)CC1=NCCN1",
    "C1CCC#CCC1",
    "CC(=C)C1CC=C(C(=O)C1)C",
    "OC(=O)C1(C)CCC=CC1(C)C",
    "NN=C1C=CC(=O)C(=C1)C(=O)O",
    "O=S1OCCO1",
    "O=Nn1[nH]c(c(c1=O)C)C",
    "O=C1NC(=N)N(C1)C",
    "N=c1oc(=N)oc(=N)o1",
    "CCc1nc(N)nc(c1c1ccc(cc1)Cl)N",
    "CCC1(C(=O)NC(=O)NC1=O)c1ccccc1",
    "OCCc1cncc(c1)CO.Cl",
    "O=C(N(C)C)Oc1ccc[n+](c1)Cc1ccccc1.[Br-]",
    "CN1CCN(CC1)C(c1ccc(cc1)Cl)c1ccccc1.Cl.Cl",
    "Nc1ccccc1C[n+]1cscc1C.[Cl-].Cl",
    "[O-][n+]1ccccc1",
    "C1CCC2C(C1)CCCC2",
    "c1ccc2c(c1)cccc2",
    "c1ccc2c(c1)c[nH]c2",
    "c1ccc2c(c1)occ2",
    "C1CNC2N1C=CS2.[Cl-]",
    "C1=CC=C2C(=CC=C2)C1",
    "OC(=O)c1ccc2c(c1)cccc2",
    "COc1ccc2c(c1)[nH]cc2",
    "OC(=O)c1cc2ccccc2cc1Cl",
    "CC1=c2ccccc2=CC1(C)C",
    "OC(=O)c1ccnc2c1nccc2",
    "C1CCc2c(C1)nccc2",
    "C1CCC2C(C1)N=CC=C2",
    "C1CCC2=CCN=C2C1",
    "C1CCC2=C(C1)CNCC2",
    "C1=CC=C2C(=NC=CS2)C1",
    "C1=NCc2c(S1)cccc2",
    "C1CNc2[n+]1ccs2.[Cl-]",
    "N1Cn2c(C1)cnc2",
    "C1CCC2=C(C1)SN=CN2",
    "C1CNCCc2c(CCC1)cccc2",
    "BrC1=C2C=NC=C2C=N1",
    "C1CCc2c(C1)cccc2",
    "Oc1ccc(c2c1cccc2)O",
    "ClC1CCC(c2c1cccc2)O",
    "C1=CCC2=C(C1)OCC2",
    "c1ccc2c(c1)c1ccccc1cc2",
    "c1ccc2c(c1)cc1c(c2)cccc1",
    "C1=CCc2c(=C1)ccc1=Nc3c(-c21)c1ccccc1cc3",
    "c1ccc2c(c1)c1c[nH]cc1cc2",
    "O=C1N2CCCCC2Cc2c1cc1ccccc1c2",
    "c1ccc2c(c1)cc1c(c2)cc[nH]1",
    "Clc1ccc2-c3c(C(=O)c2c1)cc(cc3)[N+](=O)[O-]",
    "c1oc2c(c1)cc1c(c2)ccc2c1cc1ccc3c(c1c2)cco3",
    "CN(CC[N+]1=C2CC=CC=C2Sc2c1cccc2)C",
    "CN(CCN1c2ccccc2Cc2c1cccc2)C",
    "Nc1ccc2c(c1)[n+](C)c1c(c2)ccc(c1)N.[Cl-].Cl",
    "c1ccc2c(c1)c1C=CC=c1c1=CC=Cc21",
    "c1ccc2c(c1)CC1N(C2)CCc2c1cccc2",
    "O=C1CCC2(C(=C1)CCC1C2CCC2(C1CCC2)C)C",
    "OCc1ccc(c2c1C(=O)c1cccc(c1C2=O)O)O",
    "OC1CCC2(C(C1(C)C)CCC1(C2CC=C2C1(C)CCC1(C2CC(C)(C)CC1)C)C)C",
    "Oc1ccc2c(c1)OCC1(C2c2cc(O)c(cc2C1)O)O",
    "c1cc2cccc3c2c(c1)C=C3",
    "c1cc2ccc3c4c2c(c1)ccc4ccc3",
    "c1ccc2c(c1)cc1c3c2C=Cc3ccc1",
    "c1cc2cc3cccc4c3c3c2c(c1)NC=c3cc4",
    "c1cc2c3c(c1)c1ccccc1c3c1c3c2c2ccccc2c3ccc1",
    "C1CC2CC3CCCCC3C3C2C(C1)c1c2c3cccc2cc2c1cccc2",
    "c1cc2ccc3c4c2c(c1)ccc4c1c2c3ccc3c2c(cc1)ccc3",
    "c1ccc2c(c1)cc1c(c2)cc2c3c1ccc1c3c(c3c2cccc3)ccc1",
    "c1cc2Cc3cccc4c3c3c2c(c1)c1nc2ccccc2nc1c3cc4",
    "C1=CC=c2c(C1)cc1c3c2c2ccccc2cc3ccc1",
    "c1cc2cc3cccc4c3c3c2c(c1)c1cccc2c1c3c(O4)cc2",
    "c1ccc2c(c1)cc1c(c2)c2c3cccc4c3c(c3c2c(c1)ccc3)ccc4",
    "c1ccc2c(c1)c1sc3c4c1c1c2cccc1c1c4c(c2c3cccc2)ccc1",
    "C1CC2C3C(C4C2C(C1)CCC4)C1C3C2CCCC3C2C1CCC3",
    "c1cc2cc3C=CCc4c3c3C2c(c1)cc1c3c(c4)ccc1",
    "C1CC2CC3CCCC4C3C35C2C(C1)CC1C3C(CCC1)CC1C5C(C4)CCC1",
    "C1CC2CC3CCCC4C3C35C2C(C1)CC1C3C(C2C5C(C4)CCC2)CCC1",
    "C1=CC2=CC=CC3=CN4C(=Nc5cccc6c5c4ccc6)C(=C1)C23",
    "C1=CC2=CC=CC3=c4c(C(=C1)C23)c1c(=c2c3=C1CC=Cc3ccc2)s4",
    "c1cc2c3ccc4c5c3c(c3c2c2c1N=Cc2cc3)ccc5N=C4",
    "C1CC2CC3CCCC4C3C3C2C(C1)C1CC2CC5CC6CCCC7C6C6C5C(C2CC1C3CC4)CC1C6C(C7)CCC1",
    "COc1cc2c(cc1O)ncc(c2O)c1cncc(c1)C(=O)O",
    "CN1CC(O)CC(C1)c1cnccn1",
    "CN1CCCCC1c1ccccn1",
    "C1CCC(CC1)c1cc(cc2c1cccc2)c1ccc2c(c1)cccc2",
    "C1CCC(CC1)C1CCCC(C1)c1c(ccc2c1ccnc2)c1ccc2c(c1)ccnc2",
    "s1ccc(c1)c1cccs1",
    "CN(CCN(c1ccccn1)Cc1ccc(s1)Br)C",
    "CCN(CCOc1cccc(n1)C(c1ccco1)(c1ccc(cc1)c1ccccc1)O)CC",
    "n1ccc(cc1)c1cc2ccccc2cc1C1=CCOC=C1",
    "n1ccc(cc1)c1cc2ccccc2cc1C1CCOCC1",
    "Brc1c(C2CCOCC2)c(c2ccncc2)c(c2c1cccc2)Cl",
    "C1CC[Ge]2(CC1)CCOCC2",
    "OC(=O)C1=COC=C2C1=CC1(C=COC=C1)C=C2",
    "C1CCC2(C1)SCC1(CS2)CCCC1",
    "c1ccc2c(c1)O[As]1(O2)Oc2c(O1)cccc2",
    "c1ccc2c(c1)COC12c2ccc3c(c2Oc2c1ccc1c2cccc1)cccc3",
    "O=C1C=CC2(O1)CCC2(C)C",
    "C1CCC23C(C1)CCCC3CNCC2",
    "C1CC2CC1CC2",
    "C1CC2C(C1)COC2",
    "C1CC2CC1OC2",
    "C1CC2CCC1CC2",
    "C1CC2CCC1NC2",
    "OC1CC2CCC(C1)N2C",
    "O1CC2CCC(C1)CO2",
    "N1CC2COCC(C1)C2",
    "C1CC2CCC1SC2",
    "C1CC2OCC(C1)C2",
    "C1CC2CCC(O1)C2",
    "C1CCC2CC(C1)CNC2",
    "C1CCC2C(C1)C1CCC2CC1",
    "C1CC2CCC3(C1)C2CCC3",
    "C1CC2C(C1)C1C3C(C2C2C1CCCC2)CCCC3",
    "C1CC23CCCC(C1)(C2)C3",
    "C1CCC23C(C1)(CCC3)CCC2",
    "C1CCC23C(C1)CC(CC2)C1C3CCC1",
    "C1CC2C3C1C1C2C1C3",
    "c1ccc2c(c1)c1ccc2s1",
    "C1OCc2c(O1)c1ccc2[nH]1",
    "c1ccc2c(c1)CC1C2=C2C=CC1=C2",
    "C1=CCC23C(=C1)C=C(O2)c1c3coc1",
    "C1CC2COC3(C(C1)CCC3)C2",
    "C1CCC2C(C1)N1CCCC3C1=C2C1=NCCC3C1",
    "O1C=C2CC1=C1OCC34C1=C2C=C3O4",
    "C1CC2CCC3CC2C(C1)CO3",
    "C1OC2CCC3CCC1CC3C2",
    "C1CCC2C(C1)C1CCC2CC1",
    "C1CCC2C(C1)CC1C(C2)C2CCC1C1C2CCCC1",
    "C1CCC2C(C1)CC1C(C2)C2C3C(C1C1C2COC1)CCCC3",
    "C1C=CC23C(=C1)C=C(CC2)c1c3cccc1",
    "C1CCC23C(C1)C1CCC4C(C1(CC2)OO3)COC4",
    "O1CC2C(C1)C13CCCC3CC2CC1",
    "O1CC2C3(C1)C1C3CCC2C1",
    "C1CCC2C(C1)CC13C(C2CC3)CCCC1",
    "C1=CC2=C3CCN(C2=C1)C1C23C=CC1=C2",
    "C1CCC23C(C1)CC(CC2)OO3",
    "C1C2C3CCC(C2CC2C1C1CCC2C1)C3",
    "C1CC2CC1C1C2C2C(C1)C1CC2CC1",
    "C1C2C3C(C1C1C2C2CCC1C2)C1CC3CC1",
    "C1CC2CCC1C1C2C2CC1CC2",
    "C1CCC2C(C1)C1CC2C2CC1C1C2CCCC1",
    "C1CCC2C(C1)C1CC3CC(CC2C1)C1C3CCCC1",
    "C1CCC2C(C1)C1CC3CC(C2C1)C1C3CCCC1",
    "C1C2CC3CC1CC(C2)C3",
    "C1N2CN3CN1CN(C2)C3",
    "BrC12CC3CC(C2)CC(C1)C3",
    "C1CC2C3C1CC2CC3",
    "N1CC2CCC3N(C1)C2CC3",
    "C1CC2C3C1C2CC3",
    "C1CC2C3C1CC2OC3",
    "C1CC2OCC3C(C1)C2CC1C3O1",
    "C1CC2C(C1)C1C3C2C(C1)OCO3",
    "C1CC2CC3OC(C1)C2CC3",
    "C1CC2CC3CC(C1)C2CC3",
    "C1CC2CC3C(C1)C(C2)CO3",
    "C1CC2C1C1CC3C2C(C1)CO3",
    "C[N+]12OSC3CC(C2C2CCCCC32)NC2C1CCCC2",
    "C1CCC2C(C1)CC1C3C2OC(C1)CC3",
    "C1OC2CCC3CCC1CC3C2",
    "C1OC2C(O1)C1ON3CC1C(C2)CC3",
    "N1CCC2C(C1)CC1C32CC2C(C3)CCC1N2",
    "C1CC2OC3C45C2C(C1)C1CCN(C1C5CCC3)CC4",
    "C1COSON2C3CC(OC1)CC2CC3",
    "C1C2C3C2C2C1C32",
    "C1CC2C3C1C1CCC2C(C1)O3",
    "C1C2C3C4C1C1C2C2C(O3)C4C1C2",
    "C1CCC2C(C1)C1CC3C4C2C2C5C(C1CC3C42)CCCC5",
    "C1C2C3CC4C5C2CC1C5C3C4",
    "C1OC2C3C4N(C1)C2CC(O3)C4",
    "O1CC2C(C1)C1C3C2C2C1C2C3",
    "C1C2C3C4C1C1C2C2C3CC4C12",
    "C1CC2C3C1C1CC3C3C2C13",
    "C1C2C3C4C1C1C2C3CC41",
    "C1C2C3C4C1C1C2C(C3)C41",
    "C1C2C3CC4C2CN1C4C3",
    "C1C2C3CC1C1C3C3C2CC1C3",
    "N1CN2N(C1)C1C3C2C2C1C2C3",
    "O1CC2C3C4C(C1)C2CC3OC4",
    "C1CC2C3CC4OC2C(C1)C(C3)C4",
    "C1C2CC3C1C1C4C2C3C1CC4",
    "O1CC2C3C1CC1C2COC1C3",
    "C1CC2CC1C1C2C2C3C1C1C2C1C3",
    "C1CCC2C(C1)NC1C32CC2C(C3)C3CC1N2CC3",
    "c1ccc2c(c1)C=C1C3=CC(=CN1)N=Cc1ccc(SOS23)cc1",
    "O1CC2C(C1)C1C3C(C2C2C1CCCC2)CC1C(C3)C2C3C(C1C1C2COC1)CCCC3",
    "c1ccc2c(c1)cc1c(c2)c2ccc1nc1ccc(n2)c2c1cc1c(c2)cccc1",
    "c1ccc2c(c1)c1Cc3[nH]c(c4c3cccc4)Cc3[nH]c(Cc4[nH]c(Cc2[nH]1)c1ccccc41)c1c3ccc1",
    "C1Cc2cc3ccccc3cc2CCCc2c(C1)cccc2",
    "C1CCC2C(C1)CC1C(C2)CCCC2CC(CCC1)CCC2",
    "c1cc2ssc3cccc4c3cccc4ssc3c(c1)c2ccc3",
    "c1ccc2c(c1)CN1C2=Nc2ccc(-c3ccc1cc3)cc2",
    "c1ccc2c(c1)cc1c(c2)C2Oc3ccc(OC1c1c2cc2c(c1)cccc2)cc3",
    "c1cc2Oc3ccc(cc3)Cc3nccc4c3cc(Oc3c5c(Cc(c1)c2)nccc5ccc3)cc4",
    "C1CC2=NC1=Cc1ccc([nH]1)C=c1ccc(=CC3NC4=C2CCC4C3)[nH]1",
    "C1=CCC2=CC=CC2=[C-]1.C1CCC2C(C1)CC[CH-]2.[Fe+2]",
    "C1=CCC2=CC=CC2=[C-]1.CC1=[C-]C=CC1.[Fe+2]",
    "c1c[c-]c(cc1)C1=[C-]C=CC1.c1c[c-]c(cc1)C1=[C-]C=CC1.[Fe+2].[Fe+2]",
    "ClC1=CCC(=[C-]1)C1=[C-]C=CC1.CC1=CCC(=[C-]1)C1=[C-]C=CC1.[Fe+2].[Fe+2]",
    "C1C[NH2][Cu][NH2]1.[OH-].[OH-]",
    "CC1C=C(C)C(=O)[Au](O1)(C)C",
    "C1C[NH2][Cu][NH2]1.[OH-].[OH-]",
    "C1=Cc2cccc3c2N(C1)[Zn]1(O3)Oc2c3N1CC=Cc3ccc2",
    "Cl[Pt]12N3C=CC=CC3=C3N2C(=C2N1C=CC=C2)C=CC3.[Cl-]",
    "C1=CN2C(=C3N([Fe]2)C=CC=C3)C=C1.[O-]S(=O)(=O)[O-]",
    "[O-]C(=O)C.[Na+]",
    "[CH2-][N+]1(C)CCCCC1",
    "C[N+]1(C)CCCCC1.[Cl-]",
    "c1ccc2c(c1)[S+]CS2.[Cl-]",
    "C1CC[N+]2(CC1)Cc1c(C2)cccc1.[Br-]",
    "[O-]C(=O)C(=O)[O-].[Zn+2]",
    "C1CC[N+]2(CC1)C=CC1=CC=CC1=N2"
  };
 
  fprintf(stderr, "WLN Read Benchmark\n"); 
  
  std::string buffer; 
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  
  conv.SetOutFormat("can");
  conv.AddOption("h",OpenBabel::OBConversion::OUTOPTIONS);
  
  unsigned int n_correct = 0; 
  for (unsigned int i=0; i<BENCH_N;i++) {
    if (!C_ReadWLN(wln_bench[i], &mol))
      fprintf(stderr, "%s null read\n", smiles_bench[i]); 
    else {
      buffer = conv.WriteString(&mol,true);
      
      if (strcmp(smiles_bench[i], buffer.c_str()) != 0) {
        fprintf(stderr,"%s != %s\t%s\n",wln_bench[i], smiles_bench[i], buffer.c_str()); 
      }
      else 
        n_correct++; 
      
      buffer.clear(); 
      mol.Clear(); 
    }
  }
  
  fprintf(stderr, "%d/%d compounds correct\n", n_correct, BENCH_N); 
  fprintf(stderr, "slow benchmark completed in %d seconds\n\n", BENCH_N); 
  exit(0); 
}


static void DisplayUsage()
{
  fprintf(stderr, "readwln <options> -o<format> <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, " -h                   show the help for executable usage\n");
  fprintf(stderr, " -o                   choose output format (-osmi, -oinchi, -okey, -ocan)\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln), the parser is native C\n"
                  " with a plug in function to OpenBabel\n"
        );
  fprintf(stderr, " Input is expected to either be a WLN string\n"
                  " or a file of WLN strings seperated with newline\n"
                  " characters. Detection is done by checking for a\n"
                  " file extension seperator . \n"
        );

  DisplayUsage();
}


static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i;
  unsigned int j = 0;

  cli_inp = (const char *)0;
  format = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  for (i = 1; i < argc; i++)
  {
    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1]) {

      if (ptr[1] >= 'A' && ptr[1] <= 'Z' && !j) {
        cli_inp = ptr;
        j++; 
      }
    
      else {
        switch (ptr[1]) {

        case 'h':
          DisplayHelp();

        case 'b':
          RunBenchmark(); 

        case 'o':
          if (!strcmp(ptr, "-osmi"))
          {
            format = "smi";
            break;
          }
          else if (!strcmp(ptr, "-oinchi"))
          {
            format = "inchi";
            break;
          }
          else if(!strcmp(ptr,"-okey"))
          {
            format  = "inchikey";
            break;
          }
          else if(!strcmp(ptr,"-owln"))
          {
            format  = "WLN";
            break;
          }
          else {
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','key','wln']\n");
            DisplayUsage();
          } 
        
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
        }
      }
    }
    else {
      switch (j) {
        case 0:
          cli_inp = ptr;
          break;
        
        default:
          fprintf(stderr,"Error: wln string already set - %s\n",cli_inp);
          DisplayUsage();
      }
      j++;
    }
  }

  if (!format) {
    fprintf(stderr,"Error: no output format selected\n");
    DisplayUsage();
  }

  if (!cli_inp) {
    fprintf(stderr,"Error: no input entered\n");
    DisplayUsage();
  }

  return;
}

int main(int argc, char *argv[])
{
  FILE *fp = 0; 
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  
  ProcessCommandLine(argc, argv);
  
  conv.SetOutFormat(format);
  conv.AddOption("h",OpenBabel::OBConversion::OUTOPTIONS);
  
  if (strchr(cli_inp, '.')) {
    // treat as file 
    fp = fopen(cli_inp, "r"); 
    if (!fp) {
      fprintf(stderr,"Error: file could not be opened\n"); 
      return 1;
    }
    else {
      if (!C_ReadWLNFile(fp, &mol, &conv))
        return 1; 
      fclose(fp); 
    }
  }
  else{
    if (!C_ReadWLN(cli_inp,&mol))
      return 1;
    else
      conv.Write(&mol,&std::cout);
  }
  return 0;
}
