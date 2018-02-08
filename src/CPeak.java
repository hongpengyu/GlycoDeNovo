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
import java.util.HashSet;
import java.util.List;
import java.util.Set;
/*
Corresponding to the peak in the paper
Cpeak.InferredSuperSets is all possible candidates for this peak ---(interpret as different ion)

 */

public class CPeak implements Comparable<CPeak>{
    CSpectrum   mSpectrum;
    double      mMass = -1;  // Protonated mass
    double      mMassLow = -1;
    double      mMassHigh = -1;
    double      mIntensity = 0.000001;//normalized to sum up to 1
    double      mRawMZ = -1;
    int         mRawZ = -1;
    CPeak       mComplement;
    CPeak       mHasComplement; //HasComp means the spec has the complement originally, and complement means we added it in addComp() computationally.
    List<CTopologySuperSet> mInferredSuperSets = new ArrayList<>(); // This can be considered as the ReconstructionSet for this peak.
    List<String>            mInferredFormulas;
    List<String>            mInferredGWAFormulas;
    List<Double>            mInferredMasses;
    List<Integer>           mInferredScores;

    public CPeak(CSpectrum spectrum, double intensity, double rawMZ, int rawZ) {
        this.mSpectrum = spectrum;
        this.mIntensity = intensity;
        this.mRawMZ = rawMZ;
        this.mRawZ = rawZ;
    }

    public CPeak(CSpectrum spectrum, double mass, double intensity, double rawMZ, int rawZ) {
        this(spectrum, intensity, rawMZ, rawZ);
        this.mMass = mass;
    }

    // add mMassLow / mMassHigh
    public CPeak(CSpectrum spectrum, double mass, double intensity, double rawMZ, int rawZ, double massHigh, double massLow) {
        this(spectrum, intensity, rawMZ, rawZ);
        this.mMass = mass;
        mMassLow = massLow;
        mMassHigh = massHigh;
    }

    public CPeak(CSpectrum spectrum, double mass, double intensity, CPeak complement, double massHigh, double massLow) {
        this.mSpectrum = spectrum;
        this.mMass = mass;
        this.mIntensity = intensity;
        this.mComplement = complement;
        mMassHigh = massHigh;
        mMassLow = massLow;
    }

    /**
     * If mMass is given, compare mMass, otherwise compare mRawMZ.
     * @param peak
     * @return
     */
    @Override
    public int compareTo(CPeak peak) {
        if (peak.mMass>0 && mMass>0)
            return Double.compare(mMass, peak.mMass);
        else
            return Double.compare(mRawMZ, peak.mRawMZ);
    }

    public Set<Double> protonateRawMz() {
        double rawM = mRawMZ*mRawZ;
        double metalMass = CMass.getAtomMass(mSpectrum.mMetal);
        Set<Double> protonatedMasses = new HashSet<>();
        for (int i = 1; i <= mRawZ; i++) {
            protonatedMasses.add(rawM - i * (metalMass - CMass.Electron) - (mRawZ - i) * CMass.Proton + CMass.Proton);
        }
        return protonatedMasses;
    }

    void clearInferred() {
        clearList(mInferredSuperSets);
        clearList(mInferredFormulas);
        clearList(mInferredMasses);
        clearList(mInferredScores);
    }

    private static void clearList(List l) {
        if (l != null) {
            l.clear();
        }
    }

}
