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
// Mass of Different Atom
public class CMass {
    static final double Electron = 0.0005489;
    static final double H = 1.0078250321;
    static final double H2 = 2.0156500642;
    static final double H2O = 18.0105646863;
    static final double C = 12;
    static final double N = 14.0030740052;
    static final double O = 15.9949146221;
    static final double CH2 = 14.0156500642;
    static final double Proton = 1.007276432;
    static final double Lithium = 7.0154553836;
    static final double Sodium = 22.989769;
    static final double Cesium = 132.90545;
    static final double O18 = 2.00425;
    static final double DEUTERIUM = 17.03758;
    static final double AMINOPYRIDINE = 78.05803471;
    static final double PRAGS = 120.0687;
    public static final double permethylationMassLoss = 14.0156500642 * 2;

    public static double getAtomMass(String atom){
        switch (atom){
            case "H" : return CMass.H;
            case "Proton": return CMass.Proton;
            case "Na": return CMass.Sodium;
            case "Li": return CMass.Lithium;
            case "Cs":
            case "Cesium": return CMass.Cesium;
            case "O": return CMass.O;
            case "N": return CMass.N;
            case "C": return CMass.C;
            default: throw new InvalidParameterException("atom not found!");
        }
    }
    //
    public static double reducingEndMassCompensation(String reducingEndMod, int permethylated) {
        switch (reducingEndMod) {
            case "O18":
                return CMass.O18;
            case "Deuterium":
                return CMass.DEUTERIUM;
            case "PRAGS":
                return CMass.PRAGS;
            case "Aminopyridine":
                return CMass.AMINOPYRIDINE;
            case "Reduced":
                return CMass.H2 + permethylated * CMass.CH2;
            default:
                throw new InvalidParameterException(reducingEndMod);
        }
    }
}
