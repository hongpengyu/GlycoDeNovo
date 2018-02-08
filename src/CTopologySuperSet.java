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
public class CTopologySuperSet implements Comparable<CTopologySuperSet> {
    // Important note:
    // Candidate
    // Candidate for maybe multiple peaks not just for one single peak
    // Here TSS can represent multiple candidate for different root. Although it's mass is same,
    // it can be interpreted as different ion for different peak

    // according to the paper
    // Candidate<peakID,cmass,lmass,hmass,topoReconstructionSet,topologySet>
    /*
    peakID ---- mTargetPeaks (possible for many peaks)
    cmass --- mMassPeak : Peak mass
    lmass --- mMassLow : Lower Bound of Mass
    hmass --- mMassHigh : Upper bound of Mass
    topoReconstructionSet --- mTopologySets
    topologySet --- mTopologies
     */
    // mType : Type of root of this TSS 'B' or 'C'

    // mMassPeak : Peak mass
    // mMassLow : Lower Bound of Mass
    // mMassHigh : Upper bound of Mass
    boolean     mLegal; // true for a non empty mTopologies
    String      mType;
    double      mMassPeak = 0;
    double      mMassLow = Double.MAX_VALUE;
    double      mMassHigh = 0;
    boolean     mReconstructed = false;
    // record the peak is linked to b ion or c ion...
    Map<Integer, Integer>       mTargetPeaks = new HashMap<>(); // peak index - peak type map, 1=b, 2=c
    CGlycoDeNovo mReconstructor;
    List<CTopologySet>          mTopologySets = new ArrayList<>(); // TopologySets
    List<CTopology>             mTopologies;

    public CTopologySuperSet(String mType, double mMassPeak, CGlycoDeNovo mReconstructor, int peakIndex) {
        this.mType = mType;
        this.mMassPeak = mMassPeak;
        this.mReconstructor = mReconstructor;
        if (mType.equals("B"))
            mTargetPeaks.put(peakIndex, 1);
        else if (mType.equals("C"))
            mTargetPeaks.put(peakIndex, 2);
        else if (mType.equals("T"))
            mTargetPeaks.put(peakIndex, 3);
    }

    void addPeak(int peakIndex, int peakType) {
        if (!mTargetPeaks.containsKey(peakIndex)) {
            mTargetPeaks.put(peakIndex, peakType);
        }
    }

    void addATopoSet(CTopologySet newSet) {
        if (!newSet.mType.equals(mType))
            return;
        mTopologySets.add(newSet);
        mMassLow = Math.min(newSet.mMassLow, mMassLow);
        mMassHigh = Math.max(newSet.mMassHigh, mMassHigh);

        newSet.mTargetPeaks = new HashMap<>(mTargetPeaks);// copy constructor of hashmap
    }
    // judge that if tss is contains in this TSS
    boolean contains(CTopologySuperSet tss) {
        if (tss.mMassLow < mMassLow - 0.0000001 || tss.mMassHigh > mMassHigh + 0.0000001)
            return false;
        for (CTopologySet tssTS : tss.mTopologySets) {
            boolean notFound = true;
            for (CTopologySet mTS : mTopologySets) {
                if (mTS.equals(tssTS)) {
                    notFound = false;
                    break;
                }
            }
            if (notFound)
                return false;
        }
        return true;
    }

    @Override
    public int compareTo(CTopologySuperSet o) {
        return Double.compare(mMassPeak, o.mMassPeak);
    }

//    algorithm 2
    void reconstructFormulas() {
        // construct the formula of all TSS iteratively
        // step 1 ~ 3
        if (mReconstructed)
            return;
        // mTopologySets  == TPReconstructionSet
        mTopologies = new ArrayList<>();
        Iterator<CTopologySet> iterTS = mTopologySets.iterator();
        // step 4 loop through all element in s.topoReconstructionSet
        while (iterTS.hasNext()) {
            CTopologySet mts = iterTS.next();
            // step 5 - 18
            mts.reconstructFormulas();
            if (mts.mTopologies.isEmpty())
                iterTS.remove();
            else
                // step 19
                mTopologies.addAll(mts.mTopologies);
        }

        mReconstructed = true;
        mLegal = !mTopologies.isEmpty();

        Map<String, CTopology> formulatoTp = new HashMap<>();
        Iterator<CTopology> iterTp = mTopologies.iterator();
        // remove same formula.. only save unique formula -- Topology
        while (iterTp.hasNext()) {
            CTopology mtp = iterTp.next();
            CTopology curr = formulatoTp.putIfAbsent(mtp.mFormula, mtp);
            if (curr != null) {
                curr.mSupportPeaks.addAll(mtp.mSupportPeaks);
                // need to get unique support peaks here to avoid duplicates in results
                HashSet<CPeak> temp = new HashSet<>(curr.mSupportPeaks);
                List<CPeak> uniqueSupportPeaks = new ArrayList<>(temp);
                curr.mSupportPeaks = uniqueSupportPeaks;
                curr.mScore = (double)curr.mSupportPeaks.size();
                iterTp.remove();
            }
        }
    }
}
