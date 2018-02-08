// Copyright [2018] [Pengyu Hong at Brandeis University]
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//     http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

import java.security.InvalidParameterException;
import java.util.*;
/* This is the main part of algorithm
some important note on params correspond with algorithm in the paper

Cpeak / Cpeak.mInferredSuperSets : candidate for a certain peak
mTopologySuperSets : Candidate pool
TopologySuperSet: Candidate
TopologySet : Reconstruction (Important : it's not the same as the definition in Topology)
Topology : Topology

Documentation on different parts of this Class
CGlycoDeNovo --- constructor
1. interpretPeaks --- algorithm 1

generate all possible candidates (Cpeak.mInferredSuperSets())
including find whether C-ion can be interpreted as a former B-ion
interpretApeak -- interpret current peak (algorithm 1 step 4 - 5) here only consider 2 branch
testAndAddTopologySet --- step 6 in algorithm 1
insertIntoCurrentTSS --- step 7 in algorithm 1 -- insert all possible TSS to  candidate and pool
(Cpeak.mInferredSuperSets() and Cpeak.mInferredSuperSets()) (including B-ion and C-ion)
appendNLinkedRoot --- finally add some missing part to the root

2. reconstructFormulas --- algorithm 2
call function for final peak
1.find all unique TopologySuperSet that can be interpreted related to final peak
2.from small to large --- reconstruct final peak (step 4 - 19)
call function TSS.reconstructformular() to reconstruct all TSS

 */
public class CGlycoDeNovo {
    // All cMonosaccharideClasses
    static String[] cMonosaccharideClasses = {"Xyl", "Fuc", "Hex", "HexA", "HexNAc", "Kdo", "NeuAc", "NeuGc"};
    private CMonosaccharide[] cMonosaccharideArrays = new CMonosaccharide[cMonosaccharideClasses.length];
    // All possible bonds
    // Will be used in reconstruction part
    static boolean[][] cLegalGlycosidicBonds = {
            {true, true, true, true, true, true, true, true},
            {true, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true},
            {false, false, true, true, true, true, true, true}
    };

    double  mMassAccuracyPPM = 5; // ppm = 1000000*mass_error/exact_mass; mass_error = mMassAccuracyPPM*exact_mass/1000000;
    double  mMassAccuracyPP = 0.0000005; // ppm / 1000000
    double  mMassAccuracyDalton = 0.005; // Da
    int mMaxBranchingNum = 2; // default bi-branching
    double  mIonMass = 0; // mass of the metal Ion
    String  mIonMetal; // type of mIonMetal
    double  mIntensityThreshold = 0; // Lowest Threshold for peak Intensity
    boolean mNLinked = false; // whether is Nlinked
    int     mPermethylated = 0; // Permethylated means attach CH2 on molecular
    String  mReducingEndMod;
    double  mReducingEndCompensation = 0;
    double  mFinalPeakCompensation = 0;
    int[]   mCompositionCountThreshold = new int[8];
    List<CPeak>             mPeaks;
    List<CTopologySuperSet> mTopologySuperSets;
    boolean mCheckMinus2H = false;
    boolean mCheckGap = false;
//    double[] mPeakWeights;

    String mPossibleMonosaccharideClasses;
    // Mass Compensation for Specific experiment condition
    private Delta           mDelta;
    private Delta2          mDelta2;
    private int             mNumPeaks = 0;
    private int             mCurrentPeakIdx = -1; //in matlab it is 0, but it is java and index starts from 0
    private boolean         mFinalPeak = false;
    private double          mCurrentTargetMass = 0;
    private double          mCurrentTargetMassLow = 0;
    private double          mCurrentTargetMassHigh = 0;
    private int             mCurrentTPSuperSetSize = 0;
    private boolean         mTryCIon = false;
    private boolean         mLeafPeak = false;
    // TSSB TSSC is for no Check2H // TSSB2 TSSC2 for check2H
    private CTopologySuperSet   mCurrentTopologySuperSetB;
    private CTopologySuperSet   mCurrentTopologySuperSetB2;
    private CTopologySuperSet   mCurrentTopologySuperSetC;
    private CTopologySuperSet   mCurrentTopologySuperSetC2;
    private CTopologySuperSet[] mCurrentBranches; // TSS Branches
    // initialize
    public CGlycoDeNovo(double massAccuracyPPM) {
        mMassAccuracyPPM = massAccuracyPPM;
        mMassAccuracyPP = mMassAccuracyPPM / 1000000;
        Arrays.fill(mCompositionCountThreshold, 8000);
    }

    public CGlycoDeNovo(double massAccuracyPPM, boolean checkMinus2H, boolean checkGap) {
        mCheckMinus2H = checkMinus2H;
        mCheckGap = checkGap;
        mMassAccuracyPPM = massAccuracyPPM;
        mMassAccuracyPP = mMassAccuracyPPM / 1000000;
        Arrays.fill(mCompositionCountThreshold, 8000);
    }

    // Algorithm part 1
    // Intepret ALl Peaks
    public void interpretPeaks(CSpectrum spec) {
        // initialize all params related to spec
        if (spec.mProtonated == 1) {
            mIonMetal = "Proton";
            mIonMass = CMass.Proton;
        } else {
            mIonMetal = spec.mMetal;
            switch (spec.mMetal) {
                case "Lithium":
                    mIonMass = CMass.Lithium - CMass.Electron;

                case "Na":
                case "Sodium":
                    mIonMass = CMass.Sodium - CMass.Electron;
                    break;
                default:
                    throw new InvalidParameterException("ionMetal not found");
            }
        }
//      default params
        mMassAccuracyDalton = spec.mPrecursor * mMassAccuracyPP;
        mPermethylated = spec.mPermethylated;
        mNLinked = spec.mNLinked;
        mReducingEndMod = spec.mReducingEndMod;
        if (mReducingEndMod != null)
            mReducingEndCompensation = CMass.reducingEndMassCompensation(mReducingEndMod, mPermethylated);
        mFinalPeakCompensation = mReducingEndCompensation + mPermethylated * CMass.CH2 + CMass.H2O;
        for (int i = 0; i < cMonosaccharideClasses.length; i++) {
            cMonosaccharideArrays[i] = new CMonosaccharide(cMonosaccharideClasses[i], mPermethylated);
        }
        mDelta = new Delta(cMonosaccharideClasses.length);
        mDelta2 = new Delta2(cMonosaccharideClasses.length);
        // clear all list
        spec.clearInferred();
        mPeaks = spec.mPeaks;
        mNumPeaks = mPeaks.size();
        mCurrentPeakIdx = -1;
        // step 1 - default candidate pool
        mTopologySuperSets = new ArrayList<>();
        mFinalPeak = false;

        for (int i = 0; i < mPeaks.size(); i++) {  // step 2
            // step 3 initialize the candidate set / mass
            CPeak cPeak = mPeaks.get(i);
            if (cPeak.mIntensity < mIntensityThreshold)
                continue;
            mCurrentTargetMass = cPeak.mMass;
            mCurrentTargetMassLow = cPeak.mMassLow;
            mCurrentTargetMassHigh = cPeak.mMassHigh;

            mCurrentTopologySuperSetB = null;
            mCurrentTopologySuperSetB2 = null;
            mCurrentTopologySuperSetC = null;
            mCurrentTopologySuperSetC2 = null;
            // Judge whether the mass of the peak is too small...
            // It should be larger than a single monosaccharide
            if (mPermethylated == 1) {
                if (mCurrentTargetMass < 175)
                    continue;
            } else if (mCurrentTargetMass < 131)
                    continue;

            mCurrentPeakIdx = i;
            mCurrentTPSuperSetSize = mTopologySuperSets.size();
            // initialize branches
            mCurrentBranches = new CTopologySuperSet[4];

            // Test if the current peak can be interpreted as a C ion
            // If this is the last peak, we should set mTryCIon and mFinalPeak...
            // In this condition, we should not try C Ion and we should tell the interpretApeak(), addCurrentTSSToPool() this is the last peak
            if (mCurrentPeakIdx == mNumPeaks - 1) {
                mTryCIon = false;
                mFinalPeak = true;
            } else {
                mTryCIon = true;
                // Bion + 1 H20 == Cion
                double bIonMass = mCurrentTargetMass - CMass.H2O;
                for (int j = mCurrentTPSuperSetSize - 1; j >= 0; j--) {
                    CTopologySuperSet ctss = mTopologySuperSets.get(j);
                    // Skip if ctss is of a C ion
                    double avgMass = (ctss.mMassLow + ctss.mMassHigh) / 2;
                    // search for avgMass in BionMass - 0.01 ~ BionMass + 0.01
                    if (bIonMass > avgMass + 0.01) // 0.01 is the mass accuracy, should be changed to PPM based.
                        break;
                    else if (bIonMass >= avgMass - 0.01) {
                        // Interpret as C ion of previous candidate
                        ctss.addPeak(mCurrentPeakIdx, 2);
                        mPeaks.get(i).mInferredSuperSets.add(ctss);
                        mTryCIon = false;
                        break;
                    }
                }
            }
//          step 4 ~ 5
            interpretAPeak();

            if (i == mNumPeaks - 1 && mCurrentTopologySuperSetB == null && mCurrentTopologySuperSetC == null)
                appendNLinkedRoot(); // for final peak, add Nlinked part
            addCurrentTSSToPool(); // add CurrentTSS to Pool
        }
        System.out.println("CGlycoDeNovo::reconstruct_topology done!");
    }

    // add NLinkedRoot part for Final Peak
    private void appendNLinkedRoot() {
        double mass_error = mMassAccuracyPP * mPeaks.get(mPeaks.size() - 1).mMass;
        double targetMLow = mPeaks.get(mPeaks.size() - 1).mMass - mass_error;
        double targetMHigh = mPeaks.get(mPeaks.size() - 1).mMass + mass_error;
        double dMassBIonNoFuc = CMonosaccharideSet.GlcNAc.sacPermethylated.mass * 2
                - (CMass.CH2 * 2 * mPermethylated + CMass.H2O) * 2
                + mFinalPeakCompensation;
//        % check if there is a ion corresponding to Fuc
        CTopologySuperSet ssFuc = null;
        for (CTopologySuperSet mTSS : mTopologySuperSets) {
            for (CTopologySet cs : mTSS.mTopologySets) {
                if (cs.mRootMono.mClassID == 2 && cs.mSources[0] == null) {
                    ssFuc = mTSS;
                    break;
                }
            }
            if (ssFuc != null)
                break;
        }
        double dMassBIonWithFuc = 0;
        if (ssFuc != null) {
            dMassBIonWithFuc = dMassBIonNoFuc
                    + CMonosaccharideSet.dHex.sacPermethylated.mass
                    - (CMass.CH2 * 2 * mPermethylated + CMass.H2O);
        }

        for (int cidx = mTopologySuperSets.size() - 1; cidx >= 0; cidx--) {
            CTopologySuperSet mTSS = mTopologySuperSets.get(cidx);
            double tMassLow = mTSS.mMassLow + dMassBIonNoFuc;
            double tMassHigh = mTSS.mMassHigh + dMassBIonNoFuc;

            if (tMassHigh < targetMLow - 2 && ssFuc == null) // % too light, no need to continue.
                break;
            if (targetMHigh > tMassLow && targetMLow < tMassHigh) {
                insertIntoCurrentTSS(new CTopologySet("T", mCurrentTargetMass, tMassLow, tMassHigh, mDelta.B2B.unit[4], mDelta.B2B.unit[4], new CTopologySuperSet[]{mTSS, null, null, null}, this));
            }

            if (ssFuc != null) {
                tMassLow = mTSS.mMassLow + dMassBIonWithFuc;
                tMassHigh = mTSS.mMassHigh + dMassBIonWithFuc;

                if (targetMHigh > tMassLow && targetMLow < tMassHigh) {
                    insertIntoCurrentTSS(new CTopologySet("T", mCurrentTargetMass, tMassLow, tMassHigh,
                            mDelta.B2B.unit[4], mDelta.B2B.unit[4], new CTopologySuperSet[]{ssFuc, mTSS, null, null}, this));
                    insertIntoCurrentTSS(new CTopologySet("T", mCurrentTargetMass, tMassLow, tMassHigh,
                            mDelta.B2B.unit[4], mDelta.B2B.unit[4], 1, new CTopologySuperSet[]{ssFuc, mTSS, null, null}, this));
                }
            }
        }
    }
// add Current TSSB TSSC to Pool
    private void addCurrentTSSToPool() {
//      % add mCurrentTopologySuperSetB/B1/B2/C/C1/C2 to the search space for the peaks followed

        CPeak curPeak = mPeaks.get(mCurrentPeakIdx);
        // add C2 to TSS if not null
        if (mCurrentTopologySuperSetC2 != null) {
            if (mTopologySuperSets.isEmpty()
                    || mCurrentTopologySuperSetC2.mMassLow > mTopologySuperSets.get(mTopologySuperSets.size() - 1).mMassLow) {
                // add to end if pool empty or TSS mass is largest
                mTopologySuperSets.add(mCurrentTopologySuperSetC2);
                curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetC2);
            } else {
                // search for its correct position if mass is not largest
                for (int i = mTopologySuperSets.size() - 1; i >= 0; i--) {
                    if (mTopologySuperSets.get(i).contains(mCurrentTopologySuperSetC2)) {
                        curPeak.mInferredSuperSets.add(mTopologySuperSets.get(i));
                        // type 22 is to indicate that this is a C-2H ion
                        mTopologySuperSets.get(i).addPeak(mCurrentPeakIdx, 22);
                        break;
                    } else if (mCurrentTopologySuperSetC2.mMassLow > mTopologySuperSets.get(i).mMassLow) {
                        mTopologySuperSets.add(i + 1, mCurrentTopologySuperSetC2);
                        curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetC2);
                        break;
                    }
                }
            }
        }
        if (mCurrentTopologySuperSetC != null) {
            // add to end if pool empty or TSS mass is largest
            if (mTopologySuperSets.isEmpty()
                    || mCurrentTopologySuperSetC.mMassLow > mTopologySuperSets.get(mTopologySuperSets.size() - 1).mMassLow) {
                mTopologySuperSets.add(mCurrentTopologySuperSetC);
                curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetC);
            } else {
                // search for its correct position if mass is not largest
                for (int i = mTopologySuperSets.size() - 1; i >= 0; i--) {
                    if (mTopologySuperSets.get(i).contains(mCurrentTopologySuperSetC)) {
                        curPeak.mInferredSuperSets.add(mTopologySuperSets.get(i));
                        // type 2 means this is a C ion
                        mTopologySuperSets.get(i).addPeak(mCurrentPeakIdx, 2);
                        break;
                    } else if (mCurrentTopologySuperSetC.mMassLow > mTopologySuperSets.get(i).mMassLow) {
                        mTopologySuperSets.add(i + 1, mCurrentTopologySuperSetC);
                        curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetC);
                        break;
                    }
                }
            }
        }

        if (mCurrentTopologySuperSetB2 != null) {
            // add to end if pool empty or TSS mass is largest

            if (mTopologySuperSets.isEmpty()
                    || mCurrentTopologySuperSetB2.mMassLow > mTopologySuperSets.get(mTopologySuperSets.size() - 1).mMassLow) {
                mTopologySuperSets.add(mCurrentTopologySuperSetB2);
                curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetB2);
            } else {
                // search for its correct position
                for (int i = mTopologySuperSets.size() - 1; i >= 0; i--) {
                    if (mTopologySuperSets.get(i).contains(mCurrentTopologySuperSetB2)) {
                        curPeak.mInferredSuperSets.add(mTopologySuperSets.get(i));
                        mTopologySuperSets.get(i).addPeak(mCurrentPeakIdx, 21);
                        break;
                    } else if (mCurrentTopologySuperSetB2.mMassLow > mTopologySuperSets.get(i).mMassLow) {
                        mTopologySuperSets.add(i + 1, mCurrentTopologySuperSetB2);
                        curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetB2);
                        break;
                    }
                }
            }
        }

        if (mCurrentTopologySuperSetB != null) {
            mTopologySuperSets.add(mCurrentTopologySuperSetB);
            curPeak.mInferredSuperSets.add(mCurrentTopologySuperSetB);
        }
    }
// step 4 - 5
    private void interpretAPeak() {
        if (mCurrentTargetMass < 438) { //interpret as a mono, hence no heavier than heaviest
            double massCompensationLow = mIonMass;
            if (mPermethylated == 1)
                massCompensationLow += CMass.CH2;

            if (mFinalPeak) { // if final, we should add some compensation
                massCompensationLow = massCompensationLow + mFinalPeakCompensation;
            }
            double massCompensationHigh = massCompensationLow * (1 + mMassAccuracyPP);
            massCompensationLow = massCompensationLow * (1 - mMassAccuracyPP);
            mLeafPeak = true;
            if (!testAndAddTopologySet(massCompensationLow, massCompensationHigh))
                return;
        }
        double[] branchMassesLow = new double[2];
        double[] branchMassesHigh = new double[2];
        double massD = mIonMass + mPermethylated * CMass.CH2;
        //each branch causes a CH2 loss to the joint monosaccharide when permethylated && obj.mDelta only considers linear structure (i.e., one branch).;
        mLeafPeak = false;
        for (int i = 0; i < mCurrentTPSuperSetSize; i++) {
            mCurrentBranches[0] = mTopologySuperSets.get(i);
            if (mFinalPeak) {
                branchMassesLow[0] = mCurrentBranches[0].mMassLow + mFinalPeakCompensation * (1 - mMassAccuracyPP);
                branchMassesHigh[0] = mCurrentBranches[0].mMassHigh + mFinalPeakCompensation * (1 + mMassAccuracyPP);
            } else {
                branchMassesLow[0] = mCurrentBranches[0].mMassLow;
                branchMassesHigh[0] = mCurrentBranches[0].mMassHigh;
            }
            if (!testAndAddTopologySet(branchMassesLow[0], branchMassesHigh[0]))
                break;
        // add first branch .....
            for (int ii = i; ii < mTopologySuperSets.size(); ii++) {
                mCurrentBranches[1] = mTopologySuperSets.get(ii);
                branchMassesLow[1] = mCurrentBranches[1].mMassLow - massD;
                branchMassesHigh[1] = mCurrentBranches[1].mMassHigh - massD;
                if (!testAndAddTopologySet(branchMassesLow[1] + branchMassesLow[0],
                        branchMassesHigh[1] + branchMassesHigh[0])) // sum of low and high
                    break;
                //todo: 3 and 4 branches cases to be added if needed
            }
            // add second branch .....
            mCurrentBranches[1] = null;
            // default second branch and go to set next first branch
            branchMassesLow[1] = 0;
            branchMassesHigh[1] = 0;
        }
        mCurrentBranches[0] = null;
    }
// algorithm 2 -- reconstruct last peak
    void reconstructFormulas() {

        CPeak lastPeak = mPeaks.get(mPeaks.size() - 1);
        HashSet<CTopologySuperSet> frontier = new HashSet<>(lastPeak.mInferredSuperSets);
        // put the last candidate's InferredSuperSets in.....
        HashSet<CTopologySuperSet> buffer = new HashSet<>(frontier);

        while (!frontier.isEmpty()) {
            ArrayList<CTopologySuperSet> newFrontier = new ArrayList<>();
            for (CTopologySuperSet tss : frontier) {
                for (CTopologySet ts : tss.mTopologySets) {
                    ts.mTargetPeaks = tss.mTargetPeaks;
                    newFrontier.addAll(Arrays.asList(ts.mSources));
                }
            }
            newFrontier.removeIf(Objects::isNull); // remove nulls
            frontier = new HashSet<>(newFrontier);
            buffer.addAll(frontier);
        }

        // find all unique TSS that relevant to last peak.
        List<CTopologySuperSet> sortedUniqueTSS = new ArrayList<>(buffer);
        Collections.sort(sortedUniqueTSS);

        // reconstruct all ts from small to large
        // begin algorithm 2 (CandidateSetReconstructor()) for last peak
        for (CTopologySuperSet ts : sortedUniqueTSS) {
            ts.reconstructFormulas();
        }


        System.out.println("CGlycoDeNovo::reconstruct_formulas done!");
        // link all valid formula of TSS to peak
        for (CPeak peak : mPeaks) {
            // get all valid peak
            peak.mInferredSuperSets.removeIf(e -> !e.mLegal);
            if (peak.mInferredSuperSets.isEmpty())
                continue;
            peak.mInferredSuperSets.forEach(e -> Collections.sort(e.mTopologies));
            // need to sort the mTopologies to make sure that all formula will be in correct order (number of peak (large -> small), String sorted order)
            Map<String, CTopology> formulaToTSS = new HashMap<>();
            // make sure there is only a unique result for each formula
            for (CTopologySuperSet superSet : peak.mInferredSuperSets) {
                for (CTopology mTopology : superSet.mTopologies) {
                    formulaToTSS.putIfAbsent(mTopology.mFormula, mTopology);
                }
            }

            peak.mInferredFormulas = new ArrayList<>(formulaToTSS.keySet());
            int size = peak.mInferredFormulas.size();
            peak.mInferredMasses = new ArrayList<>(size);
            peak.mInferredScores = new ArrayList<>(size);
            peak.mInferredGWAFormulas = new ArrayList<>(size);

            for (String formula : peak.mInferredFormulas) {
                peak.mInferredMasses.add(formulaToTSS.get(formula).mMass);
                peak.mInferredScores.add(formulaToTSS.get(formula).mSupportPeaks.size());
                peak.mInferredGWAFormulas.add("freeEnd--?" + formulaToTSS.get(formula).mGWAFormula.substring(3));
            }
        }
    }
// step 6 - 8
    boolean testAndAddTopologySet(double massCompensationLow, double massCompensationHigh) {

        // seems the influence of mass, C2C = B2B / B2C = C2B

        // % Check if the branches together are already too heavy or still too light
        double temp = mCurrentTargetMass - massCompensationLow;
        if (mPermethylated == 1) {
            if (temp < 160)
                return false;
            if (mCheckGap) {
                if (temp > 860) {
                    return true;
                }
            } else if (temp > 420)
                return true;
        } else {
            if (temp < 105)
                return false;
            if (mCheckGap) {
                if (temp > 660) {
                    return true;
                }
            } else if (temp > 335)
                return true;
        }


        // record idx and whether this is a minusH result
        List<Integer> idx = new ArrayList<>();
        List<Integer> flagMinusH = new ArrayList<>();

        // deal with B2B

        // loop through B2B Delta to find whether this peak can be interpreted as a B2B + mono
        // step 6 in the paper Algorithm 1
        // search find any correct of monosaccharide link
        for (int i = 0; i < mDelta.B2B.mass.length; i++) {
            double theoMassLow = mDelta.B2B.mass[i] + massCompensationLow;
            double theoMassHigh = mDelta.B2B.mass[i] + massCompensationHigh;
            if (mCurrentTargetMassHigh > theoMassLow && mCurrentTargetMassLow < theoMassHigh) {
                idx.add(i);
                flagMinusH.add(0);
            }
        }
        // if check2H, we search for result add 2H in current peak
        if (mCheckMinus2H) {
            for (int i = 0; i < mDelta.B2B.mass.length; i++) {
                double theoMassLow = mDelta.B2B.mass[i] + massCompensationLow;
                double theoMassHigh = mDelta.B2B.mass[i] + massCompensationHigh;
                if (mCurrentTargetMassHigh + CMass.H2 > theoMassLow && mCurrentTargetMassLow + CMass.H2 < theoMassHigh) {
                    idx.add(i);
                    flagMinusH.add(2);
                }
            }
        }
        for (int a = 0; a < idx.size(); a++) {
            int i = idx.get(a);
            CMonosaccharide newUnit = mDelta.B2B.unit[i];
            double theoMassLow = mDelta.B2B.mass[i] + massCompensationLow;
            double theoMassHigh = mDelta.B2B.mass[i] + massCompensationHigh;
            if ((newUnit.mClassID == 2 && !mLeafPeak) ||
                    mCompositionCountThreshold[newUnit.mClassID - 1] < 1)
                continue;
            String mType;
            if (mFinalPeak)
                mType = "T";
            else
                mType = "B";
            // step 7-8
            int h = flagMinusH.get(a);
            insertIntoCurrentTSS(new CTopologySet(mType,
                    mCurrentTargetMass,
                    Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h), Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h),
                    newUnit, mCurrentBranches, this, h));
        }
        boolean failed = idx.size() == 0;

        // Deal with B2C

        // search in mDelta.B2C.mass
        if (mTryCIon && !mFinalPeak) {
            idx = new ArrayList<>();
            flagMinusH = new ArrayList<>();
            for (int i = 0; i < mDelta.B2C.mass.length; i++) {
                double theoMassLow = mDelta.B2C.mass[i] + massCompensationLow;
                double theoMassHigh = mDelta.B2C.mass[i] + massCompensationHigh;
                if (mCurrentTargetMassHigh > theoMassLow && mCurrentTargetMassLow < theoMassHigh) {
                    idx.add(i);
                    flagMinusH.add(0);
                }
            }
            if (mCheckMinus2H) {
                for (int i = 0; i < mDelta.B2C.mass.length; i++) {
                    double theoMassLow = mDelta.B2C.mass[i] + massCompensationLow;
                    double theoMassHigh = mDelta.B2C.mass[i] + massCompensationHigh;
                    if (mCurrentTargetMassHigh + CMass.H2 > theoMassLow && mCurrentTargetMassLow + CMass.H2 < theoMassHigh) {
                        idx.add(i);
                        flagMinusH.add(2);
                    }
                }
            }
            for (int a = 0; a < idx.size(); a++) {
                int i = idx.get(a);
                CMonosaccharide newUnit = mDelta.B2C.unit[i];
                double theoMassLow = mDelta.B2C.mass[i] + massCompensationLow;
                double theoMassHigh = mDelta.B2C.mass[i] + massCompensationHigh;
                if ((newUnit.mClassID == 2 && !mLeafPeak) ||
                        mCompositionCountThreshold[newUnit.mClassID - 1] < 1)
                    continue;
                int h = flagMinusH.get(a);
                insertIntoCurrentTSS(new CTopologySet("C",
                        mCurrentTargetMass,
                        Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h) - CMass.H2O, Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h) - CMass.H2O,
                        newUnit, mCurrentBranches, this, h));
            }
        }
        failed = failed && (idx.size() == 0);

        // if failed and the checkGap condition is true, we checkGap to find whether there might exist a gap mono
        // this part is almost same from previous part
        if (failed && mCheckGap) {
            if (mCurrentBranches[1] == null && mCurrentBranches[0] != null) {
                boolean possible = false;
                for (CTopologySet aTPS : mCurrentBranches[0].mTopologySets) {
                    if (aTPS.mRootMono.mClassID != 2) {
                        possible = true;
                        break;
                    }
                }
                if (!possible) {
                    return true;
                }
            }
            idx = new ArrayList<>();
            flagMinusH = new ArrayList<>();
            // to be continued
//        % (2) try to extend to B-ions by adding one monosaccharide
//        % @@@@ Here. Need to deal with [massCompensationLow, massCompensationHigh]
            // loop through B2B Delta2 to find whether this peak can be interpreted as a B2B + mono
            // step 6 in the paper Algorithm 1
            for (int i = 0; i < mDelta2.B2B.len; i++) {
                double theoMassLow = mDelta2.B2B.mass[i] + massCompensationLow;
                double theoMassHigh = mDelta2.B2B.mass[i] + massCompensationHigh;
                if (mCurrentTargetMassHigh > theoMassLow && mCurrentTargetMassLow < theoMassHigh) {
                    idx.add(i);
                    flagMinusH.add(0);
                }
            }
            if (mCheckMinus2H) {
                for (int i = 0; i < mDelta2.B2B.len; i++) {
                    double theoMassLow = mDelta2.B2B.mass[i] + massCompensationLow;
                    double theoMassHigh = mDelta2.B2B.mass[i] + massCompensationHigh;
                    if (mCurrentTargetMassHigh + CMass.H2 > theoMassLow && mCurrentTargetMassLow + CMass.H2 < theoMassHigh) {
                        idx.add(i);
                        flagMinusH.add(2);
                    }
                }
            }
            for (int a = 0; a < idx.size(); a++) {
                int i = idx.get(a);
                CMonosaccharide newUnit1 = mDelta2.B2B.unit[i][0];
                CMonosaccharide newUnit2 = mDelta2.B2B.unit[i][1];
                if (newUnit1.mClassID == 2) { // Rule: Fuc can only be a branch, not in a linear substructure.
                    continue;
                }
                double theoMassLow = mDelta2.B2B.mass[i] + massCompensationLow;
                double theoMassHigh = mDelta2.B2B.mass[i] + massCompensationHigh;
                String mType;
                if (mFinalPeak)
                    mType = "T";
                else
                    mType = "B";
                // step 7-8

                int h = flagMinusH.get(a);
                CTopologySet newSet = new CTopologySet(mType,
                        mCurrentTargetMass,
                        Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h), Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h),
                        newUnit2, newUnit1, mCurrentBranches, this, h);
                insertIntoCurrentTSS(newSet);
                for (int b = 0; b < 4; b++) {
                    if (mCurrentBranches[b] == null) { // correspond to matlab code
                        break;
                    }
                    insertIntoCurrentTSS(new CTopologySet(mType,
                            mCurrentTargetMass,
                            Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h), Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h),
                            newUnit2, newUnit1, b, mCurrentBranches, this, h));
                }
            }
            // check whether it is C Ion....
            // search in mDelta.B2C.mass
            if (mTryCIon && !mFinalPeak) {
                idx = new ArrayList<>();
                flagMinusH = new ArrayList<>();
                for (int i = 0; i < mDelta2.B2C.len; i++) {
                    double theoMassLow = mDelta2.B2C.mass[i] + massCompensationLow;
                    double theoMassHigh = mDelta2.B2C.mass[i] + massCompensationHigh;
                    if (mCurrentTargetMassHigh > theoMassLow && mCurrentTargetMassLow < theoMassHigh) {
                        idx.add(i);
                        flagMinusH.add(0);
                    }
                }
                if (mCheckMinus2H) {
                    for (int i = 0; i < mDelta2.B2C.len; i++) {
                        double theoMassLow = mDelta2.B2C.mass[i] + massCompensationLow;
                        double theoMassHigh = mDelta2.B2C.mass[i] + massCompensationHigh;
                        if (mCurrentTargetMassHigh - CMass.H2 > theoMassLow && mCurrentTargetMassLow - CMass.H2 < theoMassHigh) {
                            idx.add(i);
                            flagMinusH.add(2);
                        }
                    }
                }

                for (int a = 0; a < idx.size(); a++) {
                    int i = idx.get(a);
                    CMonosaccharide newUnit1 = mDelta2.B2C.unit[i][0];
                    CMonosaccharide newUnit2 = mDelta2.B2C.unit[i][1];
                    if (newUnit1.mClassID == 2) { // Rule: Fuc can only be a branch, not in a linear substructure.
                        continue;
                    }
                    double theoMassLow = mDelta2.B2C.mass[i] + massCompensationLow;
                    double theoMassHigh = mDelta2.B2C.mass[i] + massCompensationHigh;
                    String mType;

                    int h = flagMinusH.get(a);
                    CTopologySet newSet = new CTopologySet("C",
                            mCurrentTargetMass,
                            Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h) - CMass.H2O, Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h) - CMass.H2O,
                            newUnit2, newUnit1, mCurrentBranches, this, h);
                    insertIntoCurrentTSS(newSet);

                    for (int b = 0; b < 4; b++) {
                        if (mCurrentBranches[b] == null) {
                            break;
                        }
                        insertIntoCurrentTSS(new CTopologySet("C",
                                mCurrentTargetMass,
                                Math.max(theoMassLow, mCurrentTargetMassLow + CMass.H * h) - CMass.H2O, Math.min(theoMassHigh, mCurrentTargetMassHigh + CMass.H * h) - CMass.H2O,
                                newUnit2, newUnit1, b, mCurrentBranches, this, h));
                    }
                }
            }
        }

        return true;
    }
// step 11
    private void insertIntoCurrentTSS(CTopologySet newSet) {
//        % Add a CTopologySet newSet to a CTopologySuperSet of the same mass and same type
//        % in obj.mTopologySuperSets. Create a new CTopologySuperSet if necessary.

        CTopologySuperSet tss = null;
        if (newSet.mType.equals("B") || newSet.mType.equals("T")) {
            // insert as B ions
            if (newSet.mMinusH == 0) {
                if (mCurrentTopologySuperSetB == null) {
                    mCurrentTopologySuperSetB = new CTopologySuperSet(newSet.mType, mCurrentTargetMass, this, mCurrentPeakIdx);
                }
                tss = mCurrentTopologySuperSetB;
            } else if (newSet.mMinusH == 2) {
                // let ion type be 21 23 to indicate this ion is an ion with -2H
                if (mCurrentTopologySuperSetB2 == null) {
                    mCurrentTopologySuperSetB2 = new CTopologySuperSet(newSet.mType, mCurrentTargetMass, this, mCurrentPeakIdx);
                    mCurrentTopologySuperSetB2.mTargetPeaks.put(mCurrentPeakIdx, 21);
                    if (newSet.mType.equals("T")) {
                        mCurrentTopologySuperSetB2.mTargetPeaks.put(mCurrentPeakIdx, 23);
                    }
                }
                tss = mCurrentTopologySuperSetB2;
            }
        } else if (newSet.mType.equals("C")) {
//            insert as C ion
            if (newSet.mMinusH == 0) {
                if (mCurrentTopologySuperSetC == null) {
                    mCurrentTopologySuperSetC = new CTopologySuperSet(newSet.mType, mCurrentTargetMass, this, mCurrentPeakIdx);
                }
                tss = mCurrentTopologySuperSetC;
            } else if (newSet.mMinusH == 2) {
                // let ion type be 22 to indicate -2H
                if (mCurrentTopologySuperSetC2 == null) {
                    mCurrentTopologySuperSetC2 = new CTopologySuperSet(newSet.mType, mCurrentTargetMass, this, mCurrentPeakIdx);
                    mCurrentTopologySuperSetC2.mTargetPeaks.put(mCurrentPeakIdx, 22);
                }
                tss = mCurrentTopologySuperSetC2;
            }
        } else
            throw new InvalidParameterException("wrong mtype");

        tss.addATopoSet(newSet);
    }

    // possible mass difference between two type of ion with gap
    class Delta2 {
        Link2 B2B;
        Link2 B2C;
        Link2 C2C;
        Link2 C2B;
        // numClass is cMonosaccharideClasses.length here
        private Delta2(int numClass) {
            double perMassLoss = 2 * mPermethylated * CMass.CH2; // If no the root unit, lose two CH2s when permethylated
            double linkageMassLoss = CMass.H2O + perMassLoss;
            B2B = new Link2(numClass * numClass);
            B2C = new Link2(numClass * numClass);
            C2C = new Link2(numClass * numClass);
            C2B = new Link2(numClass * numClass);
            int index = 0;
            for (int i = 0; i < numClass; i++) {
                CMonosaccharide mono1 = cMonosaccharideArrays[i];
                for (int j = 0; j < numClass; j++) {
                    CMonosaccharide mono2 = cMonosaccharideArrays[j];
                    if (!cLegalGlycosidicBonds[i][j]) {
                        continue;
                    }
                    // add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a B-ion
                    B2B.mass[index] = mono1.mMass + mono2.mMass - 2 * linkageMassLoss;
                    B2B.unit[index][0] = mono1; //mono.copy;
                    B2B.unit[index][1] = mono2;

                    // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a C-ion
                    B2C.mass[index] = mono1.mMass + mono2.mMass - linkageMassLoss + CMass.H2O;
                    B2C.unit[index][0] = mono1; // mono.copy;
                    B2C.unit[index][1] = mono2;

                    // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a B-ion
                    C2B.mass[index] = mono1.mMass + mono2.mMass - CMass.H2O - 2 * linkageMassLoss;
                    C2B.unit[index][0] = mono1; // mono.copy;
                    C2B.unit[index][1] = mono2; // mono.copy;

                    // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a C-ion
                    C2C.mass[index] = mono1.mMass - linkageMassLoss * 2;
                    C2C.unit[index][0] = mono1; // mono.copy;
                    C2C.unit[index][1] = mono2; // mono.copy;
                    index++;
                }
                B2B.len = index;
                B2C.len = index;
                C2B.len = index;
                C2C.len = index;
            }
        }
    }
    // possible mass difference between two type of ion with no gap

    class Delta {
        Link B2B;
        Link B2C;
        Link C2C;
        Link C2B;
        // numClass is cMonosaccharideClasses.length here
        private Delta(int numClass) {
            double perMassLoss = 2 * mPermethylated * CMass.CH2; // If no the root unit, lose two CH2s when permethylated
            double linkageMassLoss = CMass.H2O + perMassLoss;
            B2B = new Link(numClass);
            B2C = new Link(numClass);
            C2C = new Link(numClass);
            C2B = new Link(numClass);
            for (int i = 0; i < numClass; i++) {
                CMonosaccharide mono = cMonosaccharideArrays[i];

                // add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a B-ion
                B2B.mass[i] = mono.mMass - linkageMassLoss;
                B2B.unit[i] = mono; //mono.copy;

                // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a C-ion
                B2C.mass[i] = mono.mMass - perMassLoss;
                B2C.unit[i] = mono; // mono.copy;

                // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a B-ion
                C2B.mass[i] = mono.mMass - CMass.H2O - linkageMassLoss;
                C2B.unit[i] = mono; // mono.copy;

                // add 1 CGlycoDeNovo.cMonosaccharideClasses[k] to form a C-ion
                C2C.mass[i] = mono.mMass - CMass.H2O - perMassLoss;
                C2C.unit[i] = mono; // mono.copy;
            }
        }
    }
    // mass / Monosaccharide info of ion link for B C ion
    class Link {
        double[] mass;
        CMonosaccharide[] unit;

        private Link(int numClass) {
            mass = new double[numClass];
            unit = new CMonosaccharide[numClass];
        }
    }
    class Link2 {
        double[] mass;
        CMonosaccharide[][] unit;
        int len;
        private Link2(int numClass) {
            mass = new double[numClass];
            unit = new CMonosaccharide[numClass][2];
        }
    }
}
