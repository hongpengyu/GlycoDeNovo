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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
// Topology <mass, representation,supports>
/*
mass --- mMass
representation --- mFormula
supports --- mSupportPeaks
 */
public class CTopology implements Comparable<CTopology> {
    String      mType;
    double      mMass;
    String      mFormula;
    String      mGWAFormula;
    double      mScore; // Now it is the size of mSupportPeaks
    int         mRootMonoClassID;
    int[]       mCopositionCount = new int[8];
    int         mMinusH;
    List<CPeak> mSupportPeaks; // support peaks of this topology


    public CTopology(String formula, String gwaFormula, int ClassID, Set<CPeak> peaks, double size, double mass, String type, int[] compositionCountMerged) {
        mFormula = formula;
        mGWAFormula = gwaFormula;
        mRootMonoClassID = ClassID;
        mSupportPeaks = new ArrayList<>(peaks);
        Collections.sort(mSupportPeaks);
        mScore = size;
        mMass = mass;
        mType = type;
        mCopositionCount = compositionCountMerged;
    }


    @Override
    public int compareTo(CTopology o) {
        if (mScore != o.mScore) {
            return -Double.compare(mScore, o.mScore);
        } else {
            return mFormula.compareTo(o.mFormula);
        }
    }
}
