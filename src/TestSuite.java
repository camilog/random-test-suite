import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;
import java.util.stream.IntStream;

public class TestSuite {

    private final int n;
    private String testSequence = "";
    private double alpha = 0.05;
    private final ChiSquaredDistribution csd1 = new ChiSquaredDistribution(1);
    private final ChiSquaredDistribution csd2 = new ChiSquaredDistribution(2);
    private final NormalDistribution nd01 = new NormalDistribution();
    private int mForX3;
    private int kForX4;

    private double x1Threshold;

    private double x2Threshold;

    private double x3Threshold;

    private double x4Threshold;

    private double x5Threshold;

    private double calculateThresholdForX1() {
        return this.csd1.inverseCumulativeProbability(1 - this.alpha);
    }

    private double calculateThresholdForX2() {
        return this.csd2.inverseCumulativeProbability(1 - this.alpha);
    }

    private double calculateThresholdForX3() {
        ChiSquaredDistribution csdm = new ChiSquaredDistribution(Math.pow(2, this.mForX3) - 1);
        return csdm.inverseCumulativeProbability(1 - this.alpha);
    }

    private double calculateThresholdForX4() {
        ChiSquaredDistribution csdk = new ChiSquaredDistribution(2 * this.kForX4 - 2);
        return csdk.inverseCumulativeProbability(1 - this.alpha);
    }

    private double calculateThresholdForX5() {
        return this.nd01.inverseCumulativeProbability(1 - alpha);
    }

    public TestSuite(String sequence, int base) {
        if (base == 2) {
            this.testSequence = sequence;
        }
        if (base == 16) {
            for (int i = 0; i < sequence.length(); i++) {
                switch (sequence.charAt(i)){
                    case '0': this.testSequence = this.testSequence.concat("0000"); break;
                    case '1': this.testSequence = this.testSequence.concat("0001"); break;
                    case '2': this.testSequence = this.testSequence.concat("0010"); break;
                    case '3': this.testSequence = this.testSequence.concat("0011"); break;
                    case '4': this.testSequence = this.testSequence.concat("0100"); break;
                    case '5': this.testSequence = this.testSequence.concat("0101"); break;
                    case '6': this.testSequence = this.testSequence.concat("0110"); break;
                    case '7': this.testSequence = this.testSequence.concat("0111"); break;
                    case '8': this.testSequence = this.testSequence.concat("1000"); break;
                    case '9': this.testSequence = this.testSequence.concat("1001"); break;
                    case 'A': this.testSequence = this.testSequence.concat("1010"); break;
                    case 'B': this.testSequence = this.testSequence.concat("1011"); break;
                    case 'C': this.testSequence = this.testSequence.concat("1100"); break;
                    case 'D': this.testSequence = this.testSequence.concat("1101"); break;
                    case 'E': this.testSequence = this.testSequence.concat("1110"); break;
                    case 'F': this.testSequence = this.testSequence.concat("1111"); break;
                }
            }
        }
        this.n = this.testSequence.length();
    }

    private double frequencyTest() {
        this.x1Threshold = calculateThresholdForX1();
        int n0 = countCharInString(this.testSequence, '0');
        int n1 = countCharInString(this.testSequence, '1');
        int n = this.testSequence.length();

        // x1 should follow a Chi^2 with 1 degree of freedom if n >= 10
        // ==> Threshold value of 3.84

        return Math.pow((n0 - n1), 2) / n;

    }

    private double serialTest() {
        this.x2Threshold = calculateThresholdForX2();
        int n0 = countCharInString(this.testSequence, '0');
        int n1 = countCharInString(this.testSequence, '1');
        int n = this.testSequence.length();

        int n0Square = n0*n0;
        int n1Square = n1*n1;

        int n00 = countSequenceInString(this.testSequence, "00");
        int n01 = countSequenceInString(this.testSequence, "01");
        int n10 = countSequenceInString(this.testSequence, "10");
        int n11 = countSequenceInString(this.testSequence, "11");

        int n00Square = n00*n00;
        int n10Square = n01*n01;
        int n01Square = n10*n10;
        int n11Square = n11*n11;

        double a = ((4 * (n00Square + n01Square + n10Square + n11Square)) / (n - 1.0));
        double b = (2.0 * (n0Square + n1Square) / n);

        // x2 should follow a Chi^2 with 2 degree of freedom if n >= 21
        // ==> Threshold value of 5.99

        return a - b + 1.0;

    }

    private double pokerTest() {
        int n = this.testSequence.length();
        int m = findSuitableMForX3(n);
        this.mForX3 = m;
        this.x3Threshold = calculateThresholdForX3();
        int k = (int) Math.floor(n / m);

        String[] subSequences = divideSequenceInK(this.testSequence, k, m);

        // Extract different subSequences
        Object[] differentSubSequences = extractDifferentSubSequences(subSequences);
        int[] numberOfOccurrencesSquared = new int[differentSubSequences.length];
        for (int i = 0; i < numberOfOccurrencesSquared.length; i++) {
            numberOfOccurrencesSquared[i] = countSequenceNonOverlappingInString(this.testSequence, (String) differentSubSequences[i]);
            numberOfOccurrencesSquared[i] *= numberOfOccurrencesSquared[i];
        }

        double a = Math.pow(2, m) / k;
        double b = IntStream.of(numberOfOccurrencesSquared).sum();

        return (a * b) - k;

    }

    private double runsTest() {
        ArrayList<String> runs = getRuns(this.testSequence);

        ArrayList<String> blocks = getBlocks(runs);
        ArrayList<String> gaps = getGaps(runs);

        Dictionary<Integer, Integer> blocksSizeAndCount = getSizeAndCount(blocks);
        Dictionary<Integer, Integer> gapsSizeAndCount = getSizeAndCount(gaps);

        /*for (Enumeration<Integer> e = blocksSizeAndCount.keys(); e.hasMoreElements();) {
            int a = e.nextElement();
            System.out.println(a + "," + blocksSizeAndCount.get(a));
        }

        for (Enumeration<Integer> e = gapsSizeAndCount.keys(); e.hasMoreElements();) {
            int a = e.nextElement();
            System.out.println(a + "," + gapsSizeAndCount.get(a));
        }*/

        HashMap<Integer, Double> e;
        if (blocksSizeAndCount.size() >= gapsSizeAndCount.size()) {
            e = calculateE(this.testSequence.length(), blocksSizeAndCount);
        }
        else {
            e = calculateE(this.testSequence.length(), gapsSizeAndCount);
        }

        int k = findSuitableK(e);

        this.kForX4 = k;
        this.x4Threshold = calculateThresholdForX4();

        double firstSum = sumForX4(blocksSizeAndCount, e, k);
        double secondSum = sumForX4(gapsSizeAndCount, e, k);

        // Needs to be completed

        return firstSum + secondSum;

    }

    private double autoCorrelationTest() {
        this.x5Threshold = calculateThresholdForX5();
        int n = this.testSequence.length();
//        int d = findSuitableD(n);
        int d = 78;
        int functionA = functionA(this.testSequence, n, d);

        // x5 should follow a Normal(0, 1) if (n - d) >= 10

        return 2.0 * (functionA - ((n - d) / 2.0)) / Math.sqrt(n - d);
    }

    private double autoCorrelationTest(int d) {
        this.x5Threshold = calculateThresholdForX5();
        int n = this.testSequence.length();
        int functionA = functionA(this.testSequence, n, d);

        System.out.println("d = " + d + ", A(d) = " + functionA);

        // x5 should follow a Normal(0, 1) if (n - d) >= 10

        return 2.0 * (functionA - ((n - d) / 2.0)) / Math.sqrt(n - d);
    }

    public static void main(String[] args) {
        System.out.println("Welcome to Random Numbers Test Suite");
        String binSequence = generateRandomSequence(512);
//         String binSequence = generateSequenceWithJustOnes(512);
//         String binSequence = "1110001100010001010011101111001001001001111000110001000101001110111100100100100111100011000100010100111011110010010010011110001100010001010011101111001001001001";
//         System.out.println("Random sequence: " + binSequence);

        TestSuite testSuite = new TestSuite(binSequence, 2);
//        String hexSequence = "B67E086F8DF6016F7B695E76EDEF588B05C541E44A223768BBDADA8B4F1506A331C2E20B87C589991B4F2DE0733FCBF19313CBFE8DA14E8382B3643D199B173F";
//        TestSuite testSuite = new TestSuite(hexSequence, 16);

        double x1 = testSuite.frequencyTest();
        System.out.println("X1 = " + x1);
        if (x1 < testSuite.x1Threshold) {
            System.out.println("Frequency Test Succeeded");
        }

        double x2 = testSuite.serialTest();
        System.out.println("X2 = " + x2);
        if (x2 < testSuite.x2Threshold) {
            System.out.println("Serial Test Succeeded");
        }

        double x3 = testSuite.pokerTest();
        System.out.println("X3 = " + x3);
        if (x3 < testSuite.x3Threshold) {
            System.out.println("Poker Test Succeeded");
        }

        double x4 = testSuite.runsTest();
        System.out.println("X4 = " + x4);
        if (x4 < testSuite.x4Threshold) {
            System.out.println("Runs Test Succeeded");
        }

        double x5 = testSuite.autoCorrelationTest();
        System.out.println("X5 = " + x5);
        if (x5 < testSuite.x5Threshold && x5 > -1.0 * testSuite.x5Threshold) {
            System.out.println("AutoCorrelation Test Succeeded");
        }

    }

    private double sumForX4(Dictionary<Integer, Integer> blocksSizeAndCount, HashMap<Integer, Double> e, int k) {
        double count = 0.0;
        for (Enumeration<Integer> en = blocksSizeAndCount.keys(); en.hasMoreElements();) {
            int i = en.nextElement();
            if (i <= k) {
                count += (Math.pow(blocksSizeAndCount.get(i) - e.get(i), 2) / e.get(i));
            }
        }
        return count;
    }

    private int findSuitableK(HashMap<Integer, Double> e) {
        int k = this.testSequence.length();
        while (true) {
            if (k == 0) {
                break;
            }
            while (!e.containsKey(k)) {
                k--;
            }
            if (e.get(k) >= 5) {
                break;
            }
            else {
                k--;
            }
        }
        return k;
    }

    private HashMap<Integer, Double> calculateE(int n, Dictionary<Integer, Integer> gapsSizeAndCount) {
        HashMap<Integer, Double> e = new HashMap<>();
        for (Enumeration<Integer> enumeration = gapsSizeAndCount.keys(); enumeration.hasMoreElements();) {
            int i = enumeration.nextElement();
            e.put(i, (n - i + 3) / Math.pow(2, i+2));
        }
        return e;
    }

    private int countSequenceNonOverlappingInString(String testSequence, String differentSubSequence) {
        int k = differentSubSequence.length();
        String testSequenceAux = testSequence.substring(0, testSequence.length() - (testSequence.length() % k));
        int count = 0;
        for (int i = 0; i < testSequenceAux.length(); i = i+k) {
            boolean sequencePresent = true;
            for (int j = 0; j < k; j++) {
                if (!(testSequenceAux.charAt(j+i) == differentSubSequence.charAt(j))) {
                    sequencePresent = false;
                    break;
                }
            }
            if (sequencePresent)
                count++;
        }
        return count;
    }

    private Object[] extractDifferentSubSequences(String[] subSequences) {
        ArrayList<String> a = new ArrayList<>();
        for (String subSequence : subSequences) {
            if (!a.contains(subSequence))
                a.add(subSequence);
        }
        return a.toArray();
    }

    private Dictionary<Integer, Integer> getSizeAndCount(ArrayList<String> blocks) {
        Dictionary<Integer, Integer> a = new Hashtable<>();
        for (String block : blocks) {
            int blockLength = block.length();
            try {
                int b = a.get(blockLength);
                a.put(blockLength, b + 1);
            } catch (NullPointerException e) {
                a.put(blockLength, 1);
            }
        }
        return a;
    }

    private ArrayList<String> getGaps(ArrayList<String> runs) {
        ArrayList<String> a = new ArrayList<>();
        for (String run : runs) {
            if (run.charAt(0) == '0') {
                a.add(run);
            }
        }
        return a;
    }

    private ArrayList<String> getBlocks(ArrayList<String> runs) {
        ArrayList<String> a = new ArrayList<>();
        for (String run : runs) {
            if (run.charAt(0) == '1') {
                a.add(run);
            }
        }
        return a;
    }

    private int findSuitableD(int n) {
        // TODO: How to choose this value d?
        int d = n - 10;
        while (d > Math.floor(n/2)) {
            d--;
        }
        return d;
    }

    private int countCharInString(String a, char character) {
        int count = 0;
        for (int i = 0; i < a.length(); i++) {
            if (a.charAt(i) == character) {
                count++;
            }
        }
        return count;
    }

    private int countSequenceInString(String a, String sequence) {
        int count = 0;
        for (int i = 0; i < a.length() - 1; i++) {
            if (a.charAt(i) == sequence.charAt(0) && a.charAt(i+1) == sequence.charAt(1)) {
                count++;
            }
        }
        return count;
    }

    private int findSuitableMForX3(int n) {
        int m = n;
        while (Math.floor(n / m) < 5 * (Math.pow(2, m))) {
            m--;
        }
        return m;
    }

    private String[] divideSequenceInK(String sequence, int k, int m) {
        String[] subSequences = new String[k];

        for (int i = 0; i < subSequences.length; i++) {
            subSequences[i] = sequence.substring(m*i, m*(i+1));
        }

        return subSequences;

    }

    private ArrayList<String> getRuns(String sequence) {
        ArrayList<String> runs = new ArrayList<>();

        int startOfRun = 0;
        for (int i = 1; i < sequence.length(); i++) {
            if (sequence.charAt(i) != sequence.charAt(startOfRun)) {
                runs.add(sequence.substring(startOfRun, i));
                startOfRun = i;
            }
        }
        runs.add(sequence.substring(startOfRun));

        return runs;

    }

    private int functionA(String sequence, int n, int d) {
        int sum = 0;
        for (int i = 0; i <= n - d - 1; i++) {
            int xor = xor(sequence.charAt(i), sequence.charAt(i + d));
            sum += xor;
        }
        return sum;
    }

    private int xor(char b0, char b1) {
        if (b0 == b1) {
            return 0;
        }
        return 1;
    }

    private static String generateRandomSequence(int n) {
        String randomSequence = "";
        Random random = new Random();
        for (int i = 0; i < n; i++) {
            if (random.nextBoolean())
                randomSequence = randomSequence.concat("1");
            else
                randomSequence = randomSequence.concat("0");
        }
        return randomSequence;
    }

    private static String generateSequenceWithJustOnes(int n) {
        String randomSequence = "";
        for (int i = 0; i < n; i++) {
            randomSequence = randomSequence.concat("1");
        }
        return randomSequence;
    }

}
