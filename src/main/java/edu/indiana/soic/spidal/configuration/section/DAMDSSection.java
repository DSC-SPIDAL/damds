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
            distanceMatrixFile = getProperty(p, "DistanceMatrixFile", "distance.bin");
            weightMatrixFile = getProperty(p, "WeightMatrixFile", "weights.bin");
            labelFile = getProperty(p, "LabelFile", "labels.txt");
            initialPointsFile = getProperty(p, "InitialPointsFile", "init.txt");
            pointsFile = getProperty(p, "PointsFile", "points.txt");
            timingFile = getProperty(p, "TimingFile", "timings.txt");
            summaryFile = getProperty(p, "SummaryFile", "summary.txt");

            numberDataPoints = Integer.parseInt(getProperty(p, "NumberDataPoints", "-1"));
            targetDimension = Integer.parseInt(getProperty(p, "TargetDimension", "3"));
            distanceTransform = Double.parseDouble(getProperty(p, "DistanceTransform", "1.0"));
            threshold = Double.parseDouble(getProperty(p, "Threshold", "0.000001"));
            alpha = Double.parseDouble(getProperty(p, "Alpha", "0.95"));
            cgIter = Integer.parseInt(getProperty(p, "CGIterations", "20"));
            cgErrorThreshold = Double.parseDouble(getProperty(p, "CGErrorThreshold", "1"));
            isSammon = Boolean.parseBoolean(getProperty(p, "IsSammon", "false"));

            isBigEndian = Boolean.parseBoolean(getProperty(p, "IsBigEndian", "false"));
            isMemoryMapped = Boolean.parseBoolean(getProperty(p, "IsMemoryMapped", "true"));
        } catch (IOException e) {
            throw new RuntimeException("IO exception occurred while reading configuration properties file", e);
        }
    }

    private static String getProperty(Properties p, String name, String def) {
        String val = System.getProperty(name);
        if (val == null) {
            if (def != null) {
                val = p.getProperty(name, def);
            } else {
                val = p.getProperty(name);
            }
        }
        return val;
    }

    public String distanceMatrixFile;
    public String weightMatrixFile;
    public String labelFile;
    public String initialPointsFile;
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
                                       "Initial Points File",
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
                         initialPointsFile,
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


