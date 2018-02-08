// Copyright [2018] [Pengyu Hong at Brandeis University]
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

import java.util.*;
// Reconstruction <root,branchSet,topologySet>
/*
root --- mRootMono
branchSet --- mSources
topologySet --- mTopologies

 */
// It isn't the same as the TopologySet in the paper
//

public class CTopologySet {
    String          mType; // ion type
    double          mMassNormal;
    double          mMassLow;
    double          mMassHigh;
    CMonosaccharide mMissingMono; // the missing mono if gap exist
    CMonosaccharide mRootMono;
    CGlycoDeNovo    mReconstroctor;
    boolean         mReconstrocted;
    private int     mCandidateNum;
    boolean         mLegal; // indicate whether its legal
    int             mMinusH; // indicate whether this TS has minus2H


    // when interpret a peak there might missing a monosaccharide between peak and mroot
    // 000 - 0 - 0
    // previousPeak - missingmono - mrootmono

    Map<Integer, Integer>   mTargetPeaks = new HashMap<>(); // peak index - peak type map, 1=b, 2=c
    CTopologySuperSet[]     mSources = new CTopologySuperSet[4];// branch TSS
    List<CTopology>         mTopologies;
    List<String>            mTopologyFormulas;
    private CTopology[]     mBranchTopologies = new CTopology[4];
    CMonosaccharide[]       mGapMono = new CMonosaccharide[4]; // gap branch mono if gap exist

    // some constructors with different input
    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor) {
        mType = Type;
        mMassNormal = MassNormal;
        mMassLow = MassLow;
        mMassHigh = MassHigh;
        mRootMono = RootMono;
        mSources = Sources.clone();
        mReconstroctor = Reconstroctor;
    }

    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CMonosaccharide MissingMono, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor) {
        this(Type, MassNormal, MassLow, MassHigh, RootMono, Sources, Reconstroctor);
        mMissingMono = MissingMono;
    }

    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor, int MinusH) {
        this(Type, MassNormal, MassLow, MassHigh, RootMono, Sources, Reconstroctor);
        mMinusH = MinusH;
    }
    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CMonosaccharide MissingMono, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor, int MinusH) {
        this(Type, MassNormal, MassLow, MassHigh, RootMono, Sources, Reconstroctor);
        mMissingMono = MissingMono;
        mMinusH = MinusH;
    }

    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CMonosaccharide GapMono, int b, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor, int MinusH) {
        this(Type, MassNormal, MassLow, MassHigh, RootMono, Sources, Reconstroctor);
        mGapMono[b] = GapMono;
        mMinusH = MinusH;
    }

    public CTopologySet(String Type, double MassNormal, double MassLow, double MassHigh, CMonosaccharide RootMono, CMonosaccharide GapMono, int b, CTopologySuperSet[] Sources, CGlycoDeNovo Reconstroctor) {
        this(Type, MassNormal, MassLow, MassHigh, RootMono, Sources, Reconstroctor);
        mGapMono[b] = GapMono;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        CTopologySet that = (CTopologySet) o;
        // don't do this comparison because the massNormal is diff when the Ion type is different!!
//        if (Double.compare(that.mMassNormal, mMassNormal) != 0) return false;
        if (Double.compare(that.mMassLow, mMassLow) != 0) return false;
        if (Double.compare(that.mMassHigh, mMassHigh) != 0) return false;
        if (mRootMono != null ? !mRootMono.equals(that.mRootMono) : that.mRootMono != null) return false;

        return Arrays.equals(mSources, that.mSources) && Arrays.equals(mGapMono, that.mGapMono);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        temp = Double.doubleToLongBits(mMassNormal);
        result = (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(mMassLow);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(mMassHigh);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + (mRootMono != null ? mRootMono.hashCode() : 0);
        result = 31 * result + Arrays.hashCode(mGapMono);
        result = 31 * result + Arrays.hashCode(mSources);
        return result;
    }

    void reconstructFormulas() {
        // step 5 - 7
        if (mReconstrocted)
            return;
        //  if branches not reconstructed, reconstruct first
        // step 8 - 10
        for (int i = 0; i < mSources.length; i++) {
            CTopologySuperSet s = mSources[i];
            if (s != null) {
                if (!s.mReconstructed)
                    s.reconstructFormulas();
            } else
                break;
        }
        // reconstruct branch

        mTopologyFormulas = new ArrayList<>(2000);
        mTopologies = new ArrayList<>(2000);
        mCandidateNum = 0;
        // step 11 - 18
        recFormula(0);
        mReconstrocted = true;
        mLegal = (mCandidateNum > 0);
        if (mLegal) {
            mTopologyFormulas = mTopologyFormulas.subList(0, mCandidateNum); // java sublist toIndex is exclusive!
            mTopologies = mTopologies.subList(0, mCandidateNum);
        }
    }
    // step 11 - 18
    private void recFormula(int sourceIdx) {
        // recursively reconstruct formula
        if (sourceIdx > 3 || mSources[sourceIdx] == null) {
            // reconstruct formula when for loop to the last branch
            int branchCount = sourceIdx;
            int[] compositionCountMerged = new int[8];
            compositionCountMerged[mRootMono.mClassID - 1] = 1;
            // this peaks will be add to support peaks
            HashSet<CPeak> peaks = new HashSet<>();
            for (Integer peakIdx : mTargetPeaks.keySet()) {
                peaks.add(mReconstroctor.mPeaks.get(peakIdx));
            }
            double mass = mRootMono.mMass - CMass.H2O + CMass.Proton;
            if (mReconstroctor.mPermethylated == 1)
                mass -= CMass.CH2;

            if (mType.equals("T"))
                mass += mReconstroctor.mFinalPeakCompensation;

            String formula;
            String gwaFormula;
            // add missing mono to the formula
            if (mMissingMono == null) {
                formula = mRootMono.mClass;
                gwaFormula = "--4b1D-" + mRootMono.mClass + ",p";
            } else {
                formula = mMissingMono.mClass + " " + mRootMono.mClass;
                gwaFormula = "--4b1D-" + mRootMono.mClass + ",p--4b1D-" + mMissingMono + ",p";
                compositionCountMerged[mMissingMono.mClassID - 1] += 1;
                mass += mMissingMono.mMass - CMass.H2O;
                if (mReconstroctor.mPermethylated == 1)
                    mass -= CMass.CH2 * 2;
            }

            if (branchCount > 0) {
                // append branch formula to the root
                String[] branchFormulas = new String[branchCount];
                String[] gwaBranchFormulas = new String[branchCount];
                // judge if the bonds of mono is valid
                for (int i = 0; i < branchCount; i++) {
                    // if not valid, we should return
                    if (mMissingMono == null) {
                        if (!CGlycoDeNovo.cLegalGlycosidicBonds[mBranchTopologies[i].mRootMonoClassID - 1][mRootMono.mClassID - 1])
                            return;
                    } else {
                        if (!CGlycoDeNovo.cLegalGlycosidicBonds[mBranchTopologies[i].mRootMonoClassID - 1][mMissingMono.mClassID - 1])
                            return;
                    }
                    // if valid, we should generate formula from branch
                    mass += mBranchTopologies[i].mMass - CMass.Proton;
                    if (mReconstroctor.mPermethylated == 1) {
                        mass -= CMass.CH2;
                    }

                    peaks.addAll(mBranchTopologies[i].mSupportPeaks);
                    for (int j = 0; j < compositionCountMerged.length; j++) {
                        compositionCountMerged[j] += mBranchTopologies[i].mCopositionCount[j];
                    }
                    // generate formula if no gap exist
                    if (mGapMono[i] == null) {
                        branchFormulas[i] = mBranchTopologies[i].mFormula;
                        gwaBranchFormulas[i] = mBranchTopologies[i].mGWAFormula;
                    } else {
                        // generate formula if there is a gap
                        mass += mGapMono[i].mMass - CMass.H2O;
                        if (mReconstroctor.mPermethylated == 1) {
                            mass -= CMass.CH2 * 2;
                        }
                        branchFormulas[i] = mBranchTopologies[i].mFormula + " " + mGapMono[i].mClass;
                        gwaBranchFormulas[i] = "--4b1D-" + mGapMono[i].mClass + ",p--4b1D-" + mBranchTopologies[i].mGWAFormula + ",p";
                        compositionCountMerged[mGapMono[i].mClassID - 1] += 1;
                    }
                }

                if (mMassHigh - mMassLow < 0.00001) {
                    // Want to test if obj.mMasses(2) == obj.mMasses(3). However, due to rounding error, sometimes obj.mMasses(2) may not be exactly equal to obj.mMasses(3)
                    if (Math.abs(mass - mMassLow) > mReconstroctor.mMassAccuracyDalton)
                        return;
                } else if (mass <= mMassLow || mass >= mMassHigh)
                    return;
                if (branchCount == 1)
                    formula = branchFormulas[0] + " " + formula;
                else {
                    Arrays.sort(branchFormulas);
                    for (String s : branchFormulas) {
                        formula = "[" + s + "] " + formula;
                    }
                }
                String gwaFormulaofAllBranches = gwaBranchFormulas[0];
                for (int i = 1; i < gwaBranchFormulas.length; i++) {
                    gwaFormulaofAllBranches = "(" + gwaFormulaofAllBranches + ")" + gwaBranchFormulas[i];
                }
                gwaFormula = gwaFormula + gwaFormulaofAllBranches;
            }

            for (int i = 0; i < 8; i++) {
                if (compositionCountMerged[i] > mReconstroctor.mCompositionCountThreshold[i]) {
                    return;
                }
            }
            // avoid same formula result
            for (String mTopologyFormula : mTopologyFormulas) {
                if (mTopologyFormula.equals(formula))
                    return;
            }
            mCandidateNum++;
            // add current formula to formula list
            mTopologyFormulas.add(formula);
            // add all possible result of formula and result to mTopologies

            mTopologies.add(new CTopology(formula, gwaFormula, mRootMono.mClassID, peaks, (double)peaks.size(), mass, mType, compositionCountMerged));
        } else if (!mSources[sourceIdx].mTopologies.isEmpty()) {
            // reconstruct
            // recursion
            for (CTopology aTp : mSources[sourceIdx].mTopologies) {
                mBranchTopologies[sourceIdx] = aTp;
                recFormula(sourceIdx + 1);
            }
        }
    }
}
