package edu.indiana.soic.spidal.configuration.section;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Properties;
import java.util.stream.IntStream;

public class DAMDSSection {
    public DAMDSSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            distanceMatrixFile = p.getProperty("DistanceMatrixFile", "distance.bin");
            weightMatrixFile = p.getProperty("WeightMatrixFile", "weights.bin");
            labelFile = p.getProperty("LabelFile", "labels.txt");
            pointsFile = p.getProperty("PointsFile", "points.txt");
            timingFile = p.getProperty("TimingFile", "timings.txt");
            summaryFile = p.getProperty("SummaryFile","summary.txt");

            numberDataPoints = Integer.parseInt(p.getProperty("NumberDataPoints","-1"));
            targetDimension = Integer.parseInt(p.getProperty("TargetDimension","3"));
            distanceTransform = Double.parseDouble(p.getProperty("DistanceTransform","1.0"));
            threshold = Double.parseDouble(p.getProperty("Threshold", "0.000001"));
            alpha = Double.parseDouble(p.getProperty("Alpha", "0.95"));
            cgIter = Integer.parseInt(p.getProperty("CGIterations", "20"));
            cgErrorThreshold = Double.parseDouble(p.getProperty("CGErrorThreshold", "1"));
            isSammon = Boolean.parseBoolean(p.getProperty("IsSammon", "false"));

            isBigEndian = Boolean.parseBoolean(p.getProperty("IsBigEndian", "false"));
            isMemoryMapped = Boolean.parseBoolean(p.getProperty("IsMemoryMapped", "true"));
        } catch (IOException e) {
            throw new RuntimeException("IO exception occurred while reading configuration properties file", e);
        }
    }

    public String distanceMatrixFile;
    public String weightMatrixFile;
    public String labelFile;
    public String pointsFile;
    public String timingFile;
    public String summaryFile;

    public int numberDataPoints;
    public int targetDimension;
    public double distanceTransform;
    public double threshold;
    public double alpha;
    public int cgIter;
    public double cgErrorThreshold;
    public boolean isSammon;

    public boolean isBigEndian;
    public boolean isMemoryMapped;

    private String getPadding(int count, String prefix){
        StringBuilder sb = new StringBuilder(prefix);
        IntStream.range(0,count).forEach(i -> sb.append(" "));
        return sb.toString();
    }

    public String toString(boolean centerAligned) {
        String[] params = new String[]{"DistanceMatrixFile",
                                       "WeightMatrixFile",
                                       "Label Data File",
                                       "PointsFile",
                                       "TimingFile",
                                       "SummaryFile",
                                       "NumberDataPoints",
                                       "The Target Dimension",
                                       "Distance Transform (double)",
                                       "Threshold value",
                                       "Cooling parameter (alpha)",
                                       "CG Iterations",
                                       "CG Threshold",
                                       "Sammon mapping (boolean) ",
                                       "BigEndian (boolean)",
                                       "Memory mapped (boolean)"};
        Object[] args =
            new Object[]{distanceMatrixFile,
                         weightMatrixFile,
                         labelFile,
                         pointsFile,
                         timingFile,
                         summaryFile,
                         numberDataPoints,
                         targetDimension,
                         distanceTransform,
                         threshold, alpha,
                         cgIter,
                         cgErrorThreshold,
                         isSammon,
                         isBigEndian,
                         isMemoryMapped};

        java.util.Optional<Integer> maxLength =
            Arrays.stream(params).map(String::length).reduce(Math::max);
        if (!maxLength.isPresent()) { return ""; }
        final int max = maxLength.get();
        final String prefix = "  ";
        StringBuilder sb = new StringBuilder("Parameters...\n");
        if (centerAligned) {
            IntStream.range(0, params.length).forEach(
                i -> {
                    String param = params[i];
                    sb.append(getPadding(max - param.length(), prefix))
                      .append(param).append(": ").append(args[i]).append("\n");
                });
        }
        else {
            IntStream.range(0, params.length).forEach(
                i -> {
                    String param = params[i];
                    sb.append(prefix).append(param).append(":")
                      .append(getPadding(max - param.length(), ""))
                      .append(args[i]).append("\n");
                });
        }
        return sb.toString();
    }
}


