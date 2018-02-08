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

// Def of Monosaccharide
public class CMonosaccharide {
    int mID = -1;
    String  mSymbol = "";
    String  mName = "";
    String  mClass = "";
    int     mClassID = 0;
    String  mFormula = "";
    String  mInferredFormula = "";
    String  mConciseFormula = "";
    int     mPermethylated = 0;
    int     mReduced = 0;
    double  mMass = -1;
    int[]   mLegalLinkedInVs;
    int[]   mCarbon2Vertex;
    CMonosacVertices mVertices;

    public CMonosaccharide(String symNameClass, int permethylated) {
        for (CMonosaccharideSet mono : CMonosaccharideSet.values()) {
            if (symNameClass.equals(mono.symbol) || symNameClass.equals(mono.name) || symNameClass.equals(mono.sacClass)) {
                mSymbol = mono.symbol;
                mName = mono.name;
                mClass = mono.sacClass;
                mClassID = mono.classID;
                mLegalLinkedInVs = mono.legalLinkedInVs.clone();
                mCarbon2Vertex = mono.carbon2vertex.clone();
                mPermethylated = permethylated;
                if (mPermethylated == 1) {
                    mMass = mono.sacPermethylated.mass;
                    mFormula = mono.sacPermethylated.compostition;
                    mVertices = new CMonosacVertices(mono.sacPermethylated.vertexComposition, mono.sacPermethylated.vertexPermethylated, mono.sacPermethylated.vertexMass);
                } else {
                    mMass = mono.sacNative.mass;
                    mFormula = mono.sacNative.compostition;
                    mVertices = new CMonosacVertices(mono.sacNative.vertexComposition, mono.sacNative.vertexMass);
                }
                return;
            }
        }
        throw new InvalidParameterException(symNameClass + " is not a valid monosac!");
    }

    class CMonosacVertices {
        String[] mComposition;
        String[] mModification;
        int[] mPermethylated;
        double[] mMass;

        public CMonosacVertices(String[] composition, double[] mass) {
            this.mComposition = composition.clone();
            mMass = mass.clone();
        }

        public CMonosacVertices(String[] composition, int[] permethylated, double[] mass) {
            this(composition, mass);
            mPermethylated = permethylated.clone();
            String[] mModification = new String[mComposition.length];
        }

    }
}
