// Copyright [2018] [Pengyu Hong at Brandeis University]
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

import com.sun.org.apache.xpath.internal.operations.Bool;

import java.io.File;
// This Function is mainly for test. You should input path for data and path for result.
// It will call the CSpectrum Class to finish all the process

public class Test {
//    path for windows
//    static final String TEST_DATA_PATH = ".\\tests\\data\\";
//    static final String TEST_RES_PATH = ".\\tests\\results\\";
//    path for mac
//    static final String TEST_DATA_PATH = "./tests/data/";
//    static final String TEST_RES_PATH = "./tests/results/";
    public static void main(String[] args) {
        String TEST_DATA_PATH = ".\\tests\\data\\";
        String TEST_RES_PATH = ".\\tests\\results\\";
        boolean check2H = false;
        boolean checkGap = false;
        int massAccuracyPPM = 5;
        for (String arg : args) {
            arg = arg.trim();
            if (arg.startsWith("DATA_PATH=")) {
                TEST_DATA_PATH = arg.substring(10);
            } else if (arg.startsWith("RES_PATH=")) {
                TEST_RES_PATH = arg.substring(9);
            } else if (arg.startsWith("check2H=")) {
                check2H = Boolean.parseBoolean(arg.substring(8));
            } else if (arg.startsWith(("checkGap="))) {
                checkGap = Boolean.parseBoolean(arg.substring(9));
            } else if (arg.startsWith("PPM=")){
                massAccuracyPPM = Integer.parseInt(arg.substring(4));
            } else {
                throw new IllegalArgumentException("Your input is incorrect");
            }
        }
        System.out.println("Process begin");

        File data = new File(TEST_DATA_PATH);
        // process all file if input is a Directory
        if (data.isDirectory()) {
            for (String fileName : data.list()) {
                CSpectrum spec = CSpectrum.specProcessing(TEST_DATA_PATH + fileName);
                CGlycoDeNovo reconstructor = new CGlycoDeNovo(massAccuracyPPM, check2H, checkGap);
                reconstructor.interpretPeaks(spec);
                reconstructor.reconstructFormulas();
                spec.outputTXT(TEST_RES_PATH, check2H, checkGap);
//                spec.outputGWA(TEST_RES_PATH);
            }
        } else { // if input is a file, process it and output to test_res_path
            CSpectrum spec = CSpectrum.specProcessing(TEST_DATA_PATH);
            CGlycoDeNovo reconstructor = new CGlycoDeNovo(massAccuracyPPM, check2H, checkGap);
            reconstructor.interpretPeaks(spec);
            reconstructor.reconstructFormulas();
            spec.outputTXT(TEST_RES_PATH, check2H, checkGap);
//            spec.outputGWA(TEST_RES_PATH);
        }
        System.out.println("All Work Finished!!");
    }
}
