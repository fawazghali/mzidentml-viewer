/*
 * ProteoIDViewer.java
 *
 * Created on 15-Dec-2010, 10:34:37
 * 
 * Fawaz Ghali
 * 
 */
package uk.ac.liv.proteoidviewer;

import uk.ac.liv.proteoidviewer.util.XmlFilter;
import uk.ac.liv.proteoidviewer.util.OmssaFilter;
import uk.ac.liv.proteoidviewer.util.SourceFileFilter;
import uk.ac.liv.proteoidviewer.util.CsvFileFilter;
import com.compomics.util.gui.spectrum.DefaultSpectrumAnnotation;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableModel;
import org.jdesktop.swingx.JXTable;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Fragmentation;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.proteoidviewer.util.MzIdentMLFilter;
import uk.ac.liv.proteoidviewer.util.ProgressBarDialog;
import com.compomics.util.gui.spectrum.SpectrumPanel;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import javax.swing.*;

import javax.xml.bind.JAXBException;
import org.apache.commons.lang.StringEscapeUtils;
import uk.ac.liv.mzidlib.FalseDiscoveryRate;
import uk.ac.liv.mzidlib.MzIdentMLToCSV;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzml_wrapper.MzMlWrapper;
import uk.ac.liv.mzidconverters.Omssa2mzid;
import uk.ac.liv.mzidconverters.Tandem2mzid;
import uk.ac.liv.proteoidviewer.util.*;

/**
 *
 * @author Fawaz Ghali
 */
public class ProteoIDViewer extends javax.swing.JFrame {
    // mzIdentML file

    private MzIdentMLFilter mzIdentMLFilter = new MzIdentMLFilter();
    private OmssaFilter omssaFilter = new OmssaFilter();
    private XmlFilter xmlFilter = new XmlFilter();
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller = null;
    // GUI tables
    private JXTable proteinAmbiguityGroupTable;
    private JXTable proteinDetectionHypothesisTable;
    private JXTable spectrumIdentificationItemProteinViewTable;
    private JXTable spectrumIdentificationResultTable;
    private JXTable spectrumIdentificationItemTable;
    private JXTable peptideEvidenceTable;
    private JXTable fragmentationTable;
    private JXTable spectrumIdentificationItemTablePeptideView;
    private JXTable peptideEvidenceTablePeptideView;
    private JXTable fragmentationTablePeptideView;
    private JXTable dBSequenceTable;
    //Hashmaps to store mzIdentML data
    private Fragmentation fragmentation;
    //GUI
    private ProgressBarDialog progressBarDialog;
    private JCheckBox[] filterCheckBoxIon = new JCheckBox[12];
    private JCheckBox[] filterCheckBoxCharge = new JCheckBox[12];
    private List<String> filterListIon = new ArrayList();
    private List<String> filterListCharge = new ArrayList();
    private List<IonType> ionTypeList = null;
    private JCheckBox[] filterCheckBoxIon1 = new JCheckBox[12];
    private JCheckBox[] filterCheckBoxCharge1 = new JCheckBox[12];
    private List<String> filterListIon1 = new ArrayList();
    private List<String> filterListCharge1 = new ArrayList();
    private List<IonType> ionTypeList1 = null;
    boolean newHeadersIIProteinTable = false;
    boolean newHeadersIRTable = false;
    private SpectrumPanel spectrumPanel;
    private SpectrumPanel spectrumPanel1;
    private Vector<DefaultSpectrumAnnotation> peakAnnotation = new Vector();
    private Vector<DefaultSpectrumAnnotation> peakAnnotation1 = new Vector();
    //stats
    private FalseDiscoveryRate falseDiscoveryRate = null;
    private double falsePositiveSii;
    private double truePositiveSii;
    private double fdrSii;
    private double falsePositiveProteins;
    private double truePositiveProteins;
    private double fdrProteins;
    //Lists to store mzIdentML data
    private List<ProteinDetectionHypothesis> pDHListPassThreshold;
    private List<SpectrumIdentificationItem> spectrumIdentificationItemListForSpecificResult;
    // FDR and Stats
    private List<SpectrumIdentificationItem> sIIListPassThreshold = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListBelowThreshold = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListPassThresholdRankOne = new ArrayList();
    private List<Peptide> peptideListNonReduntant = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListIsDecoyFalse = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListIsDecoyTrue = new ArrayList();
    private boolean secondTab, thirdTab, fourthTab, fifthTab, sixthTab;
    private JMzReader jmzreader = null;
    private Map<Double, Double> peakList = new HashMap();
    private HashMap<String, String> siiSirHashMap = new HashMap();
    private HashMap<String, String> cvTermHashMap = new HashMap();
    private String fileName = "";
    private String sourceFile = "";

    /**
     * Creates new form ProteoIDViewer
     */
    public ProteoIDViewer() {
        // Look and Feel
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        // Swing init components
        initComponents();

        // protein tab
        // protein Ambiguity Group Table
        proteinAmbiguityGroupTable = new JXTable() {
        };
        proteinAmbiguityGroupTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        proteinAmbiguityGroupTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                proteinAmbiguityGroupTableMouseClicked(evt);
            }
        });
        JScrollPane jProteinAmbiguityGroupScrollPane = new JScrollPane(proteinAmbiguityGroupTable);
        jProteinAmbiguityGroupPanel.setLayout(new java.awt.BorderLayout());
        jProteinAmbiguityGroupPanel.add(jProteinAmbiguityGroupScrollPane);
        proteinAmbiguityGroupTable.getTableHeader().setReorderingAllowed(false);
        // protein Detection Hypothesis Table
        proteinDetectionHypothesisTable = new JXTable() {
        };
        proteinDetectionHypothesisTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        proteinDetectionHypothesisTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                proteinDetectionHypothesisTableMouseClicked(evt);
            }
        });

        JScrollPane jProteinDetectionHypothesisPane = new JScrollPane(proteinDetectionHypothesisTable);
        jProteinDetectionHypothesisPanel.setLayout(new java.awt.BorderLayout());
        jProteinDetectionHypothesisPanel.add(jProteinDetectionHypothesisPane);
        proteinDetectionHypothesisTable.getTableHeader().setReorderingAllowed(false);
        // spectrum Identification Item Protein View Table
        spectrumIdentificationItemProteinViewTable = new JXTable() {
        };
        spectrumIdentificationItemProteinViewTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationItemProteinViewTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationItemProteinTableeMouseClicked(evt);
            }
        });
        JScrollPane jspectrumIdentificationItemProteinTableScrollPane = new JScrollPane(spectrumIdentificationItemProteinViewTable);
        jSpectrumIdentificationItemProteinPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationItemProteinPanel.add(jspectrumIdentificationItemProteinTableScrollPane);
        spectrumIdentificationItemProteinViewTable.getTableHeader().setReorderingAllowed(false);
        // spectrum tab
        // spectrum Identification Result Table
        spectrumIdentificationResultTable = new JXTable() {
        };
        spectrumIdentificationResultTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationResultTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationResultTableMouseClicked(evt);
            }
        });
        JScrollPane jSpectrumIdentificationResultTableScrollPane = new JScrollPane(spectrumIdentificationResultTable);
        jSpectrumIdentificationResultPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationResultPanel.add(jSpectrumIdentificationResultTableScrollPane);
        spectrumIdentificationResultTable.getTableHeader().setReorderingAllowed(false);

        //spectrum Identification Item Table
        spectrumIdentificationItemTable = new JXTable() {
        };

        JScrollPane jSpectrumIdentificationItemTableScrollPane = new JScrollPane(spectrumIdentificationItemTable);
        jSpectrumIdentificationItemPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationItemPanel.add(jSpectrumIdentificationItemTableScrollPane);
        spectrumIdentificationItemTable.getTableHeader().setReorderingAllowed(false);
        spectrumIdentificationItemTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationItemTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationItemTableMouseClicked(evt);
            }
        });

        //peptide Evidence Table 
        peptideEvidenceTable = new JXTable() {
        };

        JScrollPane jPeptideEvidenceTableScrollPane = new JScrollPane(peptideEvidenceTable);
        jPeptideEvidencePanel.setLayout(new java.awt.BorderLayout());
        jPeptideEvidencePanel.add(jPeptideEvidenceTableScrollPane);
        peptideEvidenceTable.getTableHeader().setReorderingAllowed(false);
        peptideEvidenceTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        peptideEvidenceTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                peptideEvidenceTableMouseClicked(evt);
            }
        });

        //fragmentation Table
        fragmentationTable = new JXTable() {
        };
        fragmentationTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        fragmentationTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                fragmentationTableMouseClicked(evt);
            }
        });
        JScrollPane jFragmentationTableScrollPane = new JScrollPane(fragmentationTable);
        jFragmentationPanel.setLayout(new java.awt.BorderLayout());
        jFragmentationPanel.add(jFragmentationTableScrollPane);
        fragmentationTable.getTableHeader().setReorderingAllowed(false);

        // peptide view
        //spectrum Identification Item Table
        spectrumIdentificationItemTablePeptideView = new JXTable() {
        };

        // spectrumIdentificationItemTablePeptideView.setAutoCreateRowSorter(true);
        JScrollPane jSpectrumIdentificationItemTablePeptideViewScrollPane = new JScrollPane(spectrumIdentificationItemTablePeptideView);
        jSpectrumIdentificationItemPanel1.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationItemPanel1.add(jSpectrumIdentificationItemTablePeptideViewScrollPane);
        spectrumIdentificationItemTablePeptideView.getTableHeader().setReorderingAllowed(false);
        spectrumIdentificationItemTablePeptideView.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationItemTablePeptideView.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationItemTablePeptideViewMouseClicked(evt);
            }
        });

        //peptide Evidence Table 
        peptideEvidenceTablePeptideView = new JXTable() {
        };

        JScrollPane jPeptideEvidenceTableScrollPane1 = new JScrollPane(peptideEvidenceTablePeptideView);
        jPeptideEvidencePanel1.setLayout(new java.awt.BorderLayout());
        jPeptideEvidencePanel1.add(jPeptideEvidenceTableScrollPane1);
        peptideEvidenceTablePeptideView.getTableHeader().setReorderingAllowed(false);
        peptideEvidenceTablePeptideView.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        peptideEvidenceTablePeptideView.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                peptideEvidenceTablePeptideViewMouseClicked(evt);
            }

            private void peptideEvidenceTablePeptideViewMouseClicked(MouseEvent evt) {
            }
        });

        //fragmentation Table
        fragmentationTablePeptideView = new JXTable() {
        };
        fragmentationTablePeptideView.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        fragmentationTablePeptideView.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                fragmentationTablePeptideViewMouseClicked(evt);
            }

            private void fragmentationTablePeptideViewMouseClicked(MouseEvent evt) {
            }
        });
        JScrollPane jFragmentationTablePeptideViewScrollPane = new JScrollPane(fragmentationTablePeptideView);
        jFragmentationPanel1.setLayout(new java.awt.BorderLayout());
        jFragmentationPanel1.add(jFragmentationTablePeptideViewScrollPane);
        fragmentationTablePeptideView.getTableHeader().setReorderingAllowed(false);

        //dBSequenceTable Table
        dBSequenceTable = new JXTable() {
        };
        dBSequenceTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        dBSequenceTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                dBSequenceTableMouseClicked(evt);
            }
        });
        JScrollPane jdBSequenceTableScrollPane = new JScrollPane(dBSequenceTable);
        dBSequencePanel.setLayout(new java.awt.BorderLayout());
        dBSequencePanel.add(jdBSequenceTableScrollPane);
        dBSequenceTable.getTableHeader().setReorderingAllowed(false);

        // tables cannot be edited
//        proteinAmbiguityGroupTable.setEditable(false);
//        proteinDetectionHypothesisTable.setEditable(false);
//        spectrumIdentificationItemProteinViewTable.setEditable(false);
//        spectrumIdentificationResultTable.setEditable(false);
//        spectrumIdentificationItemTable.setEditable(false);
//        fragmentationTable.setEditable(false);
//        peptideEvidenceTable.setEditable(false);
//        dBSequenceTable.setEditable(false);
        spectrumIdentificationResultTable.setToolTipText("this corresponds to Spectrum Identification Result in mzIdentML");
        Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
        setSize(dim.width, dim.height - 40);
        setLocationRelativeTo(getRootPane());
        fileChooser.addChoosableFileFilter(omssaFilter);
        fileChooser.addChoosableFileFilter(xmlFilter);
        fileChooser.addChoosableFileFilter(mzIdentMLFilter);
        repaint();
    }

    private void protocolTableMouseClicked(MouseEvent evt) {
    }

    private void peptideEvidenceTableMouseClicked(MouseEvent evt) {
        int row = peptideEvidenceTable.getSelectedRow();
        if (row != -1) {
//            row = peptideEvidenceTable.convertRowIndexToModel(row);
            String db_ref = (String) peptideEvidenceTable.getValueAt(row, 6);

            int rowCount = dBSequenceTable.getModel().getRowCount();
            for (int i = 0; i < rowCount; i++) {
                if (db_ref.equals((String) dBSequenceTable.getValueAt(i, 0))) {

                    dBSequenceTable.setRowSelectionInterval(i, i);
                }

            }
        }
        if (!fourthTab) {
            loadDBSequenceTable();
            fourthTab = true;
        }
        mainTabbedPane.setSelectedIndex(3);
    }

    private void dBSequenceTableMouseClicked(MouseEvent evt) {
    }

    public void createTables() {

        jComboBox1.removeAllItems();
        // dBSequence view dBSequenceTable
        String[] dBSequenceTableHeaders = new String[]{"ID", "Accession", "Seq", "Protein Description"};
//        dBSequenceTable.setAutoCreateRowSorter(true);
        dBSequenceTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, dBSequenceTableHeaders) {
        });
        while (((DefaultTableModel) dBSequenceTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) dBSequenceTable.getModel()).removeRow(0);
        }
        // protocolTable
        String[] protocolTableHeaders = new String[]{"", ""};
        // protein view
        String[] proteinAmbiguityGroupTableHeaders = new String[]{"ID", "Name", "Protein Accessions", "Representative Protein", "Scores", "P-values", "Number of peptides", "Is Decoy", "passThreshold"};
        String[] proteinDetectionHypothesisTableHeaders = new String[]{"ID", "Accession", "Scores", "P-values", "Number of peptides", "Is Decoy", "passThreshold"};
        String[] spectrumIdentificationItemProteinViewTableHeaders = new String[]{"Peptide Sequence", "SII", "Name", "Score", "Expectation value", "passThreshold"};
        proteinAmbiguityGroupTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, proteinAmbiguityGroupTableHeaders) {
        });
//        proteinAmbiguityGroupTable.setAutoCreateRowSorter(true);
        proteinDetectionHypothesisTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, proteinDetectionHypothesisTableHeaders) {
        });
//        proteinDetectionHypothesisTable.setAutoCreateRowSorter(true);
        spectrumIdentificationItemProteinViewTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationItemProteinViewTableHeaders) {
        });
//        spectrumIdentificationItemProteinViewTable.setAutoCreateRowSorter(true);
        jProteinSequenceTextPane.setText("");
        jProteinSequenceTextPane.setText("");

        while (((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
        }

        //spectrum view
        int spectrumIdentificationResultCvParamLengs = 0;
        int spectrumIdentificationItemCvParamLengs = 0;
        String[] sir = null;
        String[] sii = null;
        List<SpectrumIdentificationResult> spectrumIdentificationResultList;

        spectrumIdentificationResultList = new ArrayList();

        Iterator<SpectrumIdentificationList> iterspectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
        while (iterspectrumIdentificationList.hasNext()) {
            SpectrumIdentificationList spectrumIdentificationList = iterspectrumIdentificationList.next();

            for (int j = 0; j < spectrumIdentificationList.getSpectrumIdentificationResult().size(); j++) {
                spectrumIdentificationResultList.add(spectrumIdentificationList.getSpectrumIdentificationResult().get(j));
            }
        }

        if (spectrumIdentificationResultList.size() > 0) {

            for (int i = 0; i < spectrumIdentificationResultList.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = spectrumIdentificationResultList.get(i);
                if (spectrumIdentificationResultCvParamLengs < spectrumIdentificationResult.getCvParam().size()) {
                    spectrumIdentificationResultCvParamLengs = spectrumIdentificationResult.getCvParam().size();
                    sir = new String[spectrumIdentificationResult.getCvParam().size()];
                    for (int j = 0; j < spectrumIdentificationResult.getCvParam().size(); j++) {
                        sir[j] = spectrumIdentificationResult.getCvParam().get(j).getName();
                    }
                }
                List<SpectrumIdentificationItem> spectrumIdentificationItemList = spectrumIdentificationResult.getSpectrumIdentificationItem();
                for (int j = 0; j < spectrumIdentificationItemList.size(); j++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemList.get(j);
                    if (spectrumIdentificationItemCvParamLengs < spectrumIdentificationItem.getCvParam().size()) {
                        spectrumIdentificationItemCvParamLengs = spectrumIdentificationItem.getCvParam().size();
                    }
                    sii = new String[spectrumIdentificationItem.getCvParam().size()];
                    for (int k = 0; k < spectrumIdentificationItem.getCvParam().size(); k++) {
                        sii[k] = spectrumIdentificationItem.getCvParam().get(k).getName();

                    }

                }
            }
        }

        String[] spectrumIdentificationResultTableHeaders = new String[spectrumIdentificationResultCvParamLengs + 2];
        spectrumIdentificationResultTableHeaders[0] = "ID";
        spectrumIdentificationResultTableHeaders[1] = "Spectrum ID";

        if (sir != null) {
            for (int i = 0; i < sir.length; i++) {
                String string = sir[i];
                spectrumIdentificationResultTableHeaders[2 + i] = string;

            }
        }

        String[] spectrumIdentificationItemTableHeaders = new String[spectrumIdentificationItemCvParamLengs + 8];
        spectrumIdentificationItemTableHeaders[0] = "ID";
        spectrumIdentificationItemTableHeaders[1] = "Peptide Sequence";
        spectrumIdentificationItemTableHeaders[2] = "Modification";
        spectrumIdentificationItemTableHeaders[3] = "Calculated MassToCharge";
        spectrumIdentificationItemTableHeaders[4] = "Experimental MassToCharge";
        spectrumIdentificationItemTableHeaders[5] = "Rank";
        spectrumIdentificationItemTableHeaders[6] = "Is Decoy";
        spectrumIdentificationItemTableHeaders[7] = "PassThreshold";
        if (sii != null) {
            for (int i = 0; i < sii.length; i++) {
                String string = sii[i];

                string = string.replaceAll("\\\\", "");

                spectrumIdentificationItemTableHeaders[8 + i] = string;
            }
        }

        //
        if (spectrumIdentificationItemTable.getColumnExt("X!Tandem:expect") != null) {
            spectrumIdentificationItemTable.getColumnExt("X!Tandem:expect").setComparator(new DoubleComparator());
        }
        if (spectrumIdentificationItemTable.getColumnExt("OMSSA:evalue") != null) {
            spectrumIdentificationItemTable.getColumnExt("OMSSA:evalue").setComparator(new DoubleComparator());
        }
        if (spectrumIdentificationItemTable.getColumnExt("mascot:expectation value") != null) {
            spectrumIdentificationItemTable.getColumnExt("mascot:expectation value").setComparator(new DoubleComparator());
        }
        if (spectrumIdentificationItemTable.getColumnExt("SEQUEST:expectation value") != null) {
            spectrumIdentificationItemTable.getColumnExt("SEQUEST:expectation value").setComparator(new DoubleComparator());
        }

        String[] peptideEvidenceTableHeaders = new String[7];
        peptideEvidenceTableHeaders[0] = "Start";
        peptideEvidenceTableHeaders[1] = "End";
        peptideEvidenceTableHeaders[2] = "Pre";
        peptideEvidenceTableHeaders[3] = "Post";
        peptideEvidenceTableHeaders[4] = "IsDecoy";
        peptideEvidenceTableHeaders[5] = "Peptide Sequence";
        peptideEvidenceTableHeaders[6] = "dBSequence_ref";

        String[] fragmentationTableHeaders = new String[]{"M/Z", "Intensity", "M Error", "Ion Type", "Charge"};

        spectrumIdentificationResultTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationResultTableHeaders) {
        });
//                spectrumIdentificationResultTable.setAutoCreateRowSorter(true);

        while (((DefaultTableModel) spectrumIdentificationResultTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).removeRow(0);
        }

        spectrumIdentificationItemTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationItemTableHeaders) {
        });
//                spectrumIdentificationItemTable.setAutoCreateRowSorter(true);

        peptideEvidenceTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, peptideEvidenceTableHeaders) {
        });
//                peptideEvidenceTable.setAutoCreateRowSorter(true);

        fragmentationTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, fragmentationTableHeaders) {
        });
//                fragmentationTable.setAutoCreateRowSorter(true);

        while (((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).removeRow(0);
        }

        while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
        }

//
        spectrumIdentificationItemTablePeptideView.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationItemTableHeaders) {
        });
        //   spectrumIdentificationItemTablePeptideView.setAutoCreateRowSorter(true);

        peptideEvidenceTablePeptideView.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, peptideEvidenceTableHeaders) {
        });
//                        peptideEvidenceTablePeptideView.setAutoCreateRowSorter(true);

        fragmentationTablePeptideView.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, fragmentationTableHeaders) {
        });
//                        peptideEvidenceTablePeptideView.setAutoCreateRowSorter(true);

        while (((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).removeRow(0);
        }

        while (((DefaultTableModel) peptideEvidenceTablePeptideView.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) peptideEvidenceTablePeptideView.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) fragmentationTablePeptideView.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) fragmentationTablePeptideView.getModel()).removeRow(0);
        }

        //graph view
        while (jGraph.getComponents().length > 0) {
            jGraph.remove(0);
        }
        while (jExperimentalFilterPanel.getComponents().length > 0) {
            jExperimentalFilterPanel.remove(0);
        }
        jGraph.validate();
        jGraph.repaint();

        //graph view
        while (jGraph1.getComponents().length > 0) {
            jGraph1.remove(0);
        }
        while (jExperimentalFilterPanel1.getComponents().length > 0) {
            jExperimentalFilterPanel1.remove(0);
        }
        jGraph.validate();
        jGraph.repaint();

        jProteinDescriptionEditorPane.setText("");

    }

    private void loadProteinAmbiguityGroupTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        String protein_accessions = "";
        pDHListPassThreshold = new ArrayList();

        Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
        List<ProteinDetectionHypothesis> proteinDetectionHypothesisList;
        while (iterProteinAmbiguityGroup.hasNext()) {
            ProteinAmbiguityGroup proteinAmbiguityGroup = iterProteinAmbiguityGroup.next();

            protein_accessions = "";
            proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
            boolean anchorProtein = false;
            String anchorProteinAccession = "";
            //mzid 1.2
            boolean leadProtein = false;
            String leadProteinAccession = "";
            boolean groupRepresentativeProtein = false;
            String groupRepresentativeProteinAccession = "";

            boolean isDecoy = false;
            String score = " ";
            String number_peptide = " ";
            boolean isPassThreshold = false;
            if (proteinDetectionHypothesisList.size() > 0) {
                for (int j = 0; j < proteinDetectionHypothesisList.size(); j++) {
                    try {
                        ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(j);

                        DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, proteinDetectionHypothesis.getDBSequenceRef());

                        if (dBSequence.getAccession() != null) {
                            protein_accessions = protein_accessions + dBSequence.getAccession() + ";";
                        }

                        if (proteinDetectionHypothesis.isPassThreshold()) {
                            pDHListPassThreshold.add(proteinDetectionHypothesis);
                        }

                        List<CvParam> cvParamList = proteinDetectionHypothesis.getCvParam();
                        for (int i = 0; i < cvParamList.size(); i++) {
                            CvParam cvParam = cvParamList.get(i);
                            if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.1.") && cvParam.getName().equals("anchor protein")) {
                                anchorProtein = true;
                                anchorProteinAccession = dBSequence.getAccession();
                                isDecoy = checkIfProteinDetectionHypothesisIsDecoy(proteinDetectionHypothesis);
                                if (proteinDetectionHypothesis.getPeptideHypothesis() != null) {
                                    number_peptide = String.valueOf(proteinDetectionHypothesis.getPeptideHypothesis().size());
                                }
                                isPassThreshold = proteinDetectionHypothesis.isPassThreshold();
                            } else if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.") && cvParam.getName().equals("group representative")) {
                                groupRepresentativeProtein = true;
                                groupRepresentativeProteinAccession = dBSequence.getAccession();
                                isDecoy = checkIfProteinDetectionHypothesisIsDecoy(proteinDetectionHypothesis);
                                if (proteinDetectionHypothesis.getPeptideHypothesis() != null) {
                                    number_peptide = String.valueOf(proteinDetectionHypothesis.getPeptideHypothesis().size());
                                }
                                isPassThreshold = proteinDetectionHypothesis.isPassThreshold();

                            } else if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.") && cvParam.getName().equals("leading protein")) {
                                leadProtein = true;
                                leadProteinAccession = dBSequence.getAccession();
                                isDecoy = checkIfProteinDetectionHypothesisIsDecoy(proteinDetectionHypothesis);
                                if (proteinDetectionHypothesis.getPeptideHypothesis() != null) {
                                    number_peptide = String.valueOf(proteinDetectionHypothesis.getPeptideHypothesis().size());
                                }
                                isPassThreshold = proteinDetectionHypothesis.isPassThreshold();

                            }
                            if (cvParam.getName().contains("score")) {
                                score = cvParam.getValue();
                            }

                        }
                    } catch (JAXBException ex) {
                        Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                    }

                }
            }
            protein_accessions = protein_accessions.substring(0, protein_accessions.length() - 1);
            if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.1.") && anchorProtein) {
                ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).addRow(new Object[]{
                    proteinAmbiguityGroup.getId(),
                    proteinAmbiguityGroup.getName(),
                    protein_accessions,
                    anchorProteinAccession,
                    roundTwoDecimals(Double.valueOf(score).doubleValue()),
                    " ",
                    Integer.valueOf(number_peptide),
                    String.valueOf(isDecoy),
                    String.valueOf(isPassThreshold)
                });
            } else if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.") && groupRepresentativeProtein) {
                ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).addRow(new Object[]{
                    proteinAmbiguityGroup.getId(),
                    proteinAmbiguityGroup.getName(),
                    protein_accessions,
                    groupRepresentativeProteinAccession,
                    roundTwoDecimals(Double.valueOf(score).doubleValue()),
                    " ",
                    Integer.valueOf(number_peptide),
                    String.valueOf(isDecoy),
                    String.valueOf(isPassThreshold)
                });

            } else if (mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.") && leadProtein) {
                ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).addRow(new Object[]{
                    proteinAmbiguityGroup.getId(),
                    proteinAmbiguityGroup.getName(),
                    protein_accessions,
                    leadProteinAccession,
                    roundTwoDecimals(Double.valueOf(score).doubleValue()),
                    " ",
                    Integer.valueOf(number_peptide),
                    String.valueOf(isDecoy),
                    String.valueOf(isPassThreshold)
                });

            } else {
                ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).addRow(new Object[]{
                    proteinAmbiguityGroup.getId(),
                    proteinAmbiguityGroup.getName(),
                    protein_accessions,
                    " ",
                    " ",
                    " ",
                    " ",
                    " ",
                    " "
                });
            }

        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void loadSpectrumIdentificationResultTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
        while (iterSpectrumIdentificationResult.hasNext()) {

            SpectrumIdentificationResult spectrumIdentificationResult = iterSpectrumIdentificationResult.next();
            ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).addRow(new String[]{
                spectrumIdentificationResult.getId(),
                spectrumIdentificationResult.getSpectrumID()
            });
            List<CvParam> cvParamListspectrumIdentificationResult = spectrumIdentificationResult.getCvParam();

            for (int s = 0; s < cvParamListspectrumIdentificationResult.size(); s++) {
                CvParam cvParam = cvParamListspectrumIdentificationResult.get(s);
                String name = cvParam.getName();
                for (int j = 0; j < spectrumIdentificationResultTable.getModel().getColumnCount(); j++) {
                    if (spectrumIdentificationResultTable.getModel().getColumnName(j).equals(name)) {
                        ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).getRowCount() - 1, j);
                    }
                }

            }

        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

    }

    private void loadPeptideTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        while (((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).removeRow(0);
        }
        siiSirHashMap.clear();

        Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
        HashMap<String, Peptide> pepMap = new HashMap();

        while (iterPeptide.hasNext()) {
            Peptide pep = iterPeptide.next();
            pepMap.put(pep.getId(), pep);

        }

        Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);

        while (iterSpectrumIdentificationResult.hasNext()) {
            try {
                SpectrumIdentificationResult spectrumIdentificationResult = iterSpectrumIdentificationResult.next();
                List<SpectrumIdentificationItem> spectrumIdentificationItemList = spectrumIdentificationResult.getSpectrumIdentificationItem();
                for (int i = 0; i < spectrumIdentificationItemList.size(); i++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemList.get(i);
                    siiSirHashMap.put(spectrumIdentificationItem.getId(), spectrumIdentificationResult.getId());

                    boolean isDecoy = checkIfSpectrumIdentificationItemIsDecoy(spectrumIdentificationItem);

                    //Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem.getPeptideRef());
                    Peptide peptide = pepMap.get(spectrumIdentificationItem.getPeptideRef());
                    //Peptide peptide = spectrumIdentificationItem.getPeptide();

                    if (peptide != null) {
                        List<Modification> modificationList = peptide.getModification();
                        Modification modification;
                        String residues = null;
                        Integer location = null;
                        String modificationName = null;
                        CvParam modificationCvParam;
                        String combine = null;
                        if (modificationList.size() > 0) {
                            modification = modificationList.get(0);
                            location = modification.getLocation();
                            if (modification.getResidues().size() > 0) {
                                residues = modification.getResidues().get(0);
                            }
                            List<CvParam> modificationCvParamList = modification.getCvParam();
                            if (modificationCvParamList.size() > 0) {
                                modificationCvParam = modificationCvParamList.get(0);
                                modificationName = modificationCvParam.getName();
                            }
                        }
                        if (modificationName != null) {
                            combine = modificationName;
                        }
                        if (residues != null) {
                            combine = combine + " on residues: " + residues;
                        }
                        if (location != null) {
                            combine = combine + " at location: " + location;
                        }
                        double calculatedMassToCharge = 0;
                        if (spectrumIdentificationItem.getCalculatedMassToCharge() != null) {
                            calculatedMassToCharge = spectrumIdentificationItem.getCalculatedMassToCharge().doubleValue();
                        }
                        int rank = 10;
                        if (psmRankValue.getSelectedIndex() == 0) {
                            rank = 1;
                        } else if (psmRankValue.getSelectedIndex() == 1) {
                            rank = 2;
                        } else if (psmRankValue.getSelectedIndex() == 2) {
                            rank = 3;
                        }

                        if (spectrumIdentificationItem.getRank() <= rank) {
                            ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).addRow(new Object[]{
                                spectrumIdentificationItem.getId(),
                                peptide.getPeptideSequence(),
                                combine,
                                roundTwoDecimals(calculatedMassToCharge),
                                roundTwoDecimals(spectrumIdentificationItem.getExperimentalMassToCharge()),
                                Integer.valueOf(spectrumIdentificationItem.getRank()),
                                isDecoy,
                                spectrumIdentificationItem.isPassThreshold()
                            });

                            List<CvParam> cvParamListspectrumIdentificationItem = spectrumIdentificationItem.getCvParam();

                            for (int s = 0; s < cvParamListspectrumIdentificationItem.size(); s++) {
                                CvParam cvParam = cvParamListspectrumIdentificationItem.get(s);
                                String accession = cvParam.getAccession();

                                if (cvParam.getName().equals("peptide unique to one protein")) {
                                    ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).setValueAt(1, ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() - 1, 8 + s);
                                } else if (accession.equals("MS:1001330")
                                        || accession.equals("MS:1001172")
                                        || accession.equals("MS:1001159")
                                        || accession.equals("MS:1001328")) {
                                    // ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).setValueAt(roundScientificNumbers(Double.valueOf(cvParam.getValue()).doubleValue()), ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() - 1, 8 + s);
                                    ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).setValueAt(roundThreeDecimals(Double.valueOf(cvParam.getValue()).doubleValue()), ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() - 1, 8 + s);
                                } else {
                                    ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationItemTablePeptideView.getModel()).getRowCount() - 1, 8 + s);
                                }
                            }
                        }
                    }
                }
            } catch (Exception ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

    }

    /**
     * Creates fragmentation Table Mouse Clicked
     */
    private void fragmentationTableMouseClicked(MouseEvent evt) {
    }

    /**
     * Creates spectrum Identification Result Table Mouse Clicked
     */
    private void spectrumIdentificationResultTableMouseClicked(MouseEvent evt) {

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationResultTable.getSelectedRow();
        if (row != -1) {
//            row = spectrumIdentificationResultTable.convertRowIndexToModel(row);
            try {
                while (((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
                }
                while (jGraph.getComponents().length > 0) {
                    jGraph.remove(0);
                }
                jGraph.validate();
                jGraph.repaint();
                spectrumIdentificationItemTable.scrollRowToVisible(0);
                String sir_id = (String) spectrumIdentificationResultTable.getModel().getValueAt(row, 0);
                SpectrumIdentificationResult spectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                if (jmzreader != null) {

                    Spectrum spectrum = null;
                    String spectrumID = spectrumIdentificationResult.getSpectrumID();
                    if (sourceFile.equals("mgf")) {
                        String spectrumIndex = spectrumID.substring(6);
                        Integer index1 = Integer.valueOf(spectrumIndex) + 1;
                        spectrum = jmzreader.getSpectrumById(index1.toString());
                    }
                    if (sourceFile.equals("mzML")) {
                        spectrum = jmzreader.getSpectrumById(spectrumID);

                    }

                    peakList = spectrum.getPeakList();

                    List<Double> mzValues;
                    if (spectrum.getPeakList() != null) {
                        mzValues = new ArrayList<Double>(spectrum.getPeakList().keySet());
                    } else {
                        mzValues = Collections.emptyList();
                    }

                    double[] mz = new double[mzValues.size()];
                    double[] intensities = new double[mzValues.size()];

                    int index = 0;
                    peakAnnotation.clear();
                    for (Double mzValue : mzValues) {
                        mz[index] = mzValue;
                        intensities[index] = spectrum.getPeakList().get(mzValue);

//                        peakAnnotation.add(
//                                new DefaultSpectrumAnnotation(
//                                mz[index],
//                                intensities[index],
//                                Color.blue,
//                                ""));
                        index++;
                    }

                    spectrumPanel = new SpectrumPanel(
                            mz,
                            intensities,
                            0.0,
                            "",
                            "");
                    spectrumPanel.setAnnotations(peakAnnotation);
                    jGraph.setLayout(new java.awt.BorderLayout());
                    jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph.add(spectrumPanel);
                    jGraph.validate();
                    jGraph.repaint();

                }

                spectrumIdentificationItemListForSpecificResult = spectrumIdentificationResult.getSpectrumIdentificationItem();
                if (spectrumIdentificationItemListForSpecificResult.size() > 0) {

                    for (int i = 0; i < spectrumIdentificationItemListForSpecificResult.size(); i++) {
                        try {
                            SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(i);
                            boolean isDecoy = checkIfSpectrumIdentificationItemIsDecoy(spectrumIdentificationItem);

                            Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem.getPeptideRef());
                            if (peptide != null) {
                                List<Modification> modificationList = peptide.getModification();
                                Modification modification;
                                String residues = null;
                                Integer location = null;
                                String modificationName = null;
                                CvParam modificationCvParam;
                                String combine = null;
                                if (modificationList.size() > 0) {
                                    modification = modificationList.get(0);
                                    location = modification.getLocation();
                                    if (modification.getResidues().size() > 0) {
                                        residues = modification.getResidues().get(0);
                                    }
                                    List<CvParam> modificationCvParamList = modification.getCvParam();
                                    if (modificationCvParamList.size() > 0) {
                                        modificationCvParam = modificationCvParamList.get(0);
                                        modificationName = modificationCvParam.getName();
                                    }
                                }
                                if (modificationName != null) {
                                    combine = modificationName;
                                }
                                if (residues != null) {
                                    combine = combine + " on residues: " + residues;
                                }
                                if (location != null) {
                                    combine = combine + " at location: " + location;
                                }
                                double calculatedMassToCharge = 0;
                                if (spectrumIdentificationItem.getCalculatedMassToCharge() != null) {
                                    calculatedMassToCharge = spectrumIdentificationItem.getCalculatedMassToCharge().doubleValue();
                                }
                                ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).addRow(new Object[]{
                                    spectrumIdentificationItem.getId(),
                                    peptide.getPeptideSequence(),
                                    combine,
                                    roundTwoDecimals(calculatedMassToCharge),
                                    roundTwoDecimals(spectrumIdentificationItem.getExperimentalMassToCharge()),
                                    Integer.valueOf(spectrumIdentificationItem.getRank()),
                                    isDecoy,
                                    spectrumIdentificationItem.isPassThreshold()
                                });

                                List<CvParam> cvParamListspectrumIdentificationItem = spectrumIdentificationItem.getCvParam();

                                for (int s = 0; s < cvParamListspectrumIdentificationItem.size(); s++) {
                                    CvParam cvParam = cvParamListspectrumIdentificationItem.get(s);
                                    String accession = cvParam.getAccession();

                                    if (cvParam.getName().equals("peptide unique to one protein")) {
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(1, ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 8 + s);
                                    } else if (accession.equals("MS:1001330")
                                            || accession.equals("MS:1001172")
                                            || accession.equals("MS:1001159")
                                            || accession.equals("MS:1001328")) {
                                        // ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(roundScientificNumbers(Double.valueOf(cvParam.getValue()).doubleValue()), ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 8 + s);
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 8 + s);
                                    } else {
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 8 + s);
                                    }
                                }

                            }
                        } catch (JAXBException ex) {
                            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                }
            } catch (JMzReaderException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private double roundTwoDecimals(double d) {
        double multipicationFactor = Math.pow(10, 2);
        return Math.round(d * multipicationFactor) / multipicationFactor;
    }

    private double roundThreeDecimals(double d) {
        double multipicationFactor = Math.pow(10, 3);
        return Math.round(d * multipicationFactor) / multipicationFactor;
    }

    private String roundScientificNumbers(double d) {
        DecimalFormat formatter = new DecimalFormat("0.#E00");
        String f = formatter.format(d);
        return f;
    }

    private boolean checkIfSpectrumIdentificationItemIsDecoy(SpectrumIdentificationItem spectrumIdentificationItem) {
        boolean result = false;

        List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
        for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
            try {
                PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
                PeptideEvidence peptiedEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                if (peptiedEvidence != null && peptiedEvidence.isIsDecoy()) {
                    result = true;
                    break;
                }
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        return result;
    }

    /**
     * Creates spectrum Identification Item Table Mouse Clicked
     */
    private void spectrumIdentificationItemTablePeptideViewMouseClicked(MouseEvent evt) {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTablePeptideView.getSelectedRow();
        if (row != -1) {
//            row = spectrumIdentificationItemTablePeptideView.convertRowIndexToModel(row);
            try {
                while (((DefaultTableModel) fragmentationTablePeptideView.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) fragmentationTablePeptideView.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) peptideEvidenceTablePeptideView.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) peptideEvidenceTablePeptideView.getModel()).removeRow(0);
                }

                fragmentationTablePeptideView.scrollRowToVisible(0);
                SpectrumIdentificationItem spectrumIdentificationItem = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationItem.class, (String) spectrumIdentificationItemTablePeptideView.getValueAt(row, 0));

//                String sir_id1 = (String)siiSirHashMap.get((String) spectrumIdentificationItemTablePeptideView.getValueAt(row, 0));
//                System.out.println("SIR");
//                System.out.println(sir_id1);
//                System.out.println("SII");
//                System.out.println((String) spectrumIdentificationItemTablePeptideView.getValueAt(row, 0));
//                System.out.println("----------------------------");
//                
                if (spectrumIdentificationItem != null) {
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                    if (peptideEvidenceRefList != null) {
                        for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                            try {
                                PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);

                                PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());

                                ((DefaultTableModel) peptideEvidenceTablePeptideView.getModel()).addRow(new Object[]{
                                    peptideEvidence.getStart(),
                                    peptideEvidence.getEnd(),
                                    peptideEvidence.getPre(),
                                    peptideEvidence.getPost(),
                                    peptideEvidence.isIsDecoy(),
                                    peptideEvidence.getPeptideRef(),
                                    peptideEvidence.getDBSequenceRef()
                                // "<html><a href=>" +peptideEvidence.getDBSequenceRef()+"</a>"
                                });
                            } catch (JAXBException ex) {
                                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                }

                fragmentation = spectrumIdentificationItem.getFragmentation();

                if (fragmentation != null) {
                    ionTypeList1 = fragmentation.getIonType();
                    if (ionTypeList1 != null) {
                        for (int i = 0; i < ionTypeList1.size(); i++) {
                            IonType ionType = ionTypeList1.get(i);
                            CvParam cvParam = ionType.getCvParam();
                            if (!filterListIon1.contains(cvParam.getName())) {
                                filterListIon1.add(cvParam.getName());
                            }
                            if (!filterListCharge1.contains(String.valueOf(ionType.getCharge()))) {
                                filterListCharge1.add(String.valueOf(ionType.getCharge()));
                            }
                            List m_mz, m_intensity, m_error;
                            m_mz = ionType.getFragmentArray().get(0).getValues();
                            m_intensity = ionType.getFragmentArray().get(1).getValues();
                            m_error = ionType.getFragmentArray().get(2).getValues();
                            String type = cvParam.getName();
                            type = type.replaceAll(" ion", "");
                            type = type.replaceAll("param: ", "");

                            if (m_mz != null && !m_mz.isEmpty()) {
                                for (int j = 0; j < m_mz.size(); j++) {
                                    ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
                                        Double.valueOf(m_mz.get(j).toString()),
                                        Double.valueOf(m_intensity.get(j).toString()),
                                        Double.valueOf(m_error.get(j).toString()),
                                        type + ionType.getIndex().get(j),
                                        Integer.valueOf(ionType.getCharge())
                                    });

                                }
                            }
                        }
                    }
//                    while (jExperimentalFilterPanel1.getComponents().length > 0) {
//                        jExperimentalFilterPanel1.remove(0);
//                    }
//                    jExperimentalFilterPanel1.validate();
//                    jExperimentalFilterPanel1.repaint();
//                    if (filterListIon1.isEmpty() && filterListCharge1.isEmpty()) {
//                        jExperimentalFilterPanel1.setLayout(new GridLayout(0, 1));
//                    } else {
//                        jExperimentalFilterPanel1.setLayout(new GridLayout(0, filterListIon1.size() + filterListCharge1.size()));
//                    }
                    // int max = 12;
//                    if (filterListIon1.size() <= 12) {
//                        max = filterListIon1.size();
//                    }
//
//                    for (int k = 0; k < max; k++) {
//                        String name = filterListIon1.get(k);
//                        filterCheckBoxIon1[k] = new JCheckBox(name);
//                        filterCheckBoxIon1[k].setSelected(true);
//                        filterCheckBoxIon1[k].addItemListener(new ItemListener() {
//
//                            public void itemStateChanged(ItemEvent E) {
//                                updateGraph1();
//                            }
//                        });
//                        String type = name;
//                        type = type.replaceFirst("frag:", "");
//                        type = type.replaceFirst("ion", "");
//                        type = type.replaceFirst("internal", "");
//
//                        filterCheckBoxIon1[k].setText(type);
//                        jExperimentalFilterPanel1.add(filterCheckBoxIon1[k]);
//                    }
//
//                    for (int k = 0; k < filterListCharge1.size(); k++) {
//                        String name = filterListCharge1.get(k);
//                        filterCheckBoxCharge1[k] = new JCheckBox(name);
//                        filterCheckBoxCharge1[k].setSelected(true);
//                        filterCheckBoxCharge1[k].addItemListener(new ItemListener() {
//
//                            public void itemStateChanged(ItemEvent E) {
//                                updateGraph1();
//                            }
//                        });
//                        String type = name;
//
//
//                        filterCheckBoxCharge1[k].setText(type);
//                       jExperimentalFilterPanel1.add(filterCheckBoxCharge1[k]);
//                    }

//                        jExperimentalFilterPanel1.repaint();
//                    jExperimentalFilterPanel1.revalidate();
                    double[] mzValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                    double[] intensityValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                    double[] m_errorValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                    peakAnnotation1.clear();
                    for (int k = 0; k < fragmentationTablePeptideView.getModel().getRowCount(); k++) {
                        mzValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 0));

                        intensityValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 1));
                        m_errorValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 2));

                        String type = (String) fragmentationTablePeptideView.getModel().getValueAt(k, 3);
                        type = type.replaceFirst("frag:", "");
                        type = type.replaceFirst("ion", "");
                        type = type.replaceFirst("internal", "");

                        peakAnnotation1.add(
                                new DefaultSpectrumAnnotation(
                                        mzValuesAsDouble[k],
                                        m_errorValuesAsDouble[k],
                                        Color.blue,
                                        type));
                    }
//                    if (mzValuesAsDouble.length > 0) {
//                        spectrumPanel1 = new SpectrumPanel(
//                                mzValuesAsDouble,
//                                intensityValuesAsDouble,
//                                spectrumIdentificationItem.getExperimentalMassToCharge(),
//                                String.valueOf(spectrumIdentificationItem.getChargeState()),
//                                spectrumIdentificationItem.getName());
//
//                        spectrumPanel1.setAnnotations(peakAnnotation1);
//                       
//
//                        while (jGraph1.getComponents().length > 0) {
//                            jGraph1.remove(0);
//                        }
//                        jGraph1.setLayout(new java.awt.BorderLayout());
//                        jGraph1.setLayout(new javax.swing.BoxLayout(jGraph1, javax.swing.BoxLayout.LINE_AXIS));
//                       
//                        jGraph1.add(spectrumPanel1);
//                        jGraph1.validate();
//                        jGraph1.repaint();
//                        this.repaint();
//                    }

                    while (jGraph1.getComponents().length > 0) {
                        jGraph1.remove(0);
                    }
                    if (jmzreader != null) {
                        try {

                            String sir_id = (String) siiSirHashMap.get((String) spectrumIdentificationItemTablePeptideView.getValueAt(row, 0));
                            SpectrumIdentificationResult spectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                            Spectrum spectrum = null;
                            String spectrumID = spectrumIdentificationResult.getSpectrumID();
                            if (sourceFile.equals("mgf")) {
                                String spectrumIndex = spectrumID.substring(6);
                                Integer index1 = Integer.valueOf(spectrumIndex) + 1;
                                spectrum = jmzreader.getSpectrumById(index1.toString());
                            }
                            if (sourceFile.equals("mzML")) {
                                spectrum = jmzreader.getSpectrumById(spectrumID);
                            }

                            peakList = spectrum.getPeakList();

                            List<Double> mzValues;
                            if (spectrum.getPeakList() != null) {
                                mzValues = new ArrayList<Double>(spectrum.getPeakList().keySet());
                            } else {
                                mzValues = Collections.emptyList();
                            }

                            double[] mz = new double[mzValues.size()];
                            double[] intensities = new double[mzValues.size()];

                            int index = 0;
                            for (Double mzValue : mzValues) {
                                mz[index] = mzValue;
                                intensities[index] = spectrum.getPeakList().get(mzValue);
                                index++;
                            }
                            spectrumPanel1 = new SpectrumPanel(
                                    mz,
                                    intensities,
                                    spectrumIdentificationItem.getExperimentalMassToCharge(),
                                    String.valueOf(spectrumIdentificationItem.getChargeState()),
                                    spectrumIdentificationItem.getName());
                            spectrumPanel1.setAnnotations(peakAnnotation1);
                            jGraph1.setLayout(new java.awt.BorderLayout());
                            jGraph1.setLayout(new javax.swing.BoxLayout(jGraph1, javax.swing.BoxLayout.LINE_AXIS));
                            jGraph1.add(spectrumPanel1);
                            jGraph1.validate();
                            jGraph1.repaint();
                        } catch (JMzReaderException ex) {
                            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                        } catch (JAXBException ex) {
                            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                        }

                    } else if (mzValuesAsDouble.length > 0) {
                        spectrumPanel1 = new SpectrumPanel(
                                mzValuesAsDouble,
                                intensityValuesAsDouble,
                                spectrumIdentificationItem.getExperimentalMassToCharge(),
                                String.valueOf(spectrumIdentificationItem.getChargeState()),
                                spectrumIdentificationItem.getName());

                        spectrumPanel1.setAnnotations(peakAnnotation1);

                        jGraph1.setLayout(new java.awt.BorderLayout());
                        jGraph1.setLayout(new javax.swing.BoxLayout(jGraph1, javax.swing.BoxLayout.LINE_AXIS));
                        jGraph1.add(spectrumPanel1);
                        jGraph1.validate();
                        jGraph1.repaint();
                        this.repaint();
                    }
                }
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    /**
     * Creates spectrum Identification Item Table Mouse Clicked
     */
    private void spectrumIdentificationItemTableMouseClicked(MouseEvent evt) {

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTable.getSelectedRow();
        if (row != -1) {
//            row = spectrumIdentificationItemTable.convertRowIndexToModel(row);
            while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
            }
            while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
                ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
            }

            fragmentationTable.scrollRowToVisible(0);
            SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(row);

            if (spectrumIdentificationItem != null) {
                List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                if (peptideEvidenceRefList != null) {
                    for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                        try {
                            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);

                            PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());

                            ((DefaultTableModel) peptideEvidenceTable.getModel()).addRow(new Object[]{
                                peptideEvidence.getStart(),
                                peptideEvidence.getEnd(),
                                peptideEvidence.getPre(),
                                peptideEvidence.getPost(),
                                peptideEvidence.isIsDecoy(),
                                peptideEvidence.getPeptideRef(),
                                peptideEvidence.getDBSequenceRef()
                            // "<html><a href=>" +peptideEvidence.getDBSequenceRef()+"</a>"
                            });
                        } catch (JAXBException ex) {
                            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                }
            }

            fragmentation = spectrumIdentificationItem.getFragmentation();

            if (fragmentation != null) {
                ionTypeList = fragmentation.getIonType();
                if (ionTypeList != null) {
                    for (int i = 0; i < ionTypeList.size(); i++) {
                        IonType ionType = ionTypeList.get(i);
                        CvParam cvParam = ionType.getCvParam();
                        if (!filterListIon.contains(cvParam.getName())) {
                            filterListIon.add(cvParam.getName());
                        }
                        if (!filterListCharge.contains(String.valueOf(ionType.getCharge()))) {
                            filterListCharge.add(String.valueOf(ionType.getCharge()));
                        }
                        List m_mz, m_intensity, m_error;
                        m_mz = ionType.getFragmentArray().get(0).getValues();
                        m_intensity = ionType.getFragmentArray().get(1).getValues();
                        m_error = ionType.getFragmentArray().get(2).getValues();
                        String type = cvParam.getName();
                        type = type.replaceAll("param: ", "");
                        type = type.replaceAll(" ion", "");

                        if (m_mz != null && !m_mz.isEmpty()) {
                            for (int j = 0; j < m_mz.size(); j++) {
                                ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                    Double.valueOf(m_mz.get(j).toString()),
                                    Double.valueOf(m_intensity.get(j).toString()),
                                    Double.valueOf(m_error.get(j).toString()),
                                    type + ionType.getIndex().get(j),
                                    Integer.valueOf(ionType.getCharge())
                                });

                            }
                        }
                    }
                }
//                while (jExperimentalFilterPanel.getComponents().length > 0) {
//                    jExperimentalFilterPanel.remove(0);
//                }
//                jExperimentalFilterPanel.validate();
//                jExperimentalFilterPanel.repaint();
//                if (filterListIon.isEmpty() && filterListCharge.isEmpty()) {
//                    jExperimentalFilterPanel.setLayout(new GridLayout(0, 1));
//                } else {
//                    jExperimentalFilterPanel.setLayout(new GridLayout(0, filterListIon.size() + filterListCharge.size()));
//                }
//                int max = 12;
//                if (filterListIon.size() <= 12) {
//                    max = filterListIon.size();
//                }
//                for (int k = 0; k < max; k++) {
//                    String name = filterListIon.get(k);
//                    filterCheckBoxIon[k] = new JCheckBox(name);
//                    filterCheckBoxIon[k].setSelected(true);
//                    filterCheckBoxIon[k].addItemListener(new ItemListener() {
//
//                        public void itemStateChanged(ItemEvent E) {
//                            updateGraph();
//                        }
//                    });
//                    String type = name;
//                    type = type.replaceFirst("frag:", "");
//                    type = type.replaceFirst("ion", "");
//                    type = type.replaceFirst("internal", "");
//
//                    filterCheckBoxIon[k].setText(type);
//                    jExperimentalFilterPanel.add(filterCheckBoxIon[k]);
//                }
//
//                for (int k = 0; k < filterListCharge.size(); k++) {
//                    String name = filterListCharge.get(k);
//                    filterCheckBoxCharge[k] = new JCheckBox(name);
//                    filterCheckBoxCharge[k].setSelected(true);
//                    filterCheckBoxCharge[k].addItemListener(new ItemListener() {
//
//                        public void itemStateChanged(ItemEvent E) {
//                            updateGraph();
//                        }
//                    });
//                    String type = name;
//
//
//                    filterCheckBoxCharge[k].setText(type);
//                    jExperimentalFilterPanel.add(filterCheckBoxCharge[k]);
//                }
//
//                jExperimentalFilterPanel.repaint();
//                jExperimentalFilterPanel.revalidate();

                double[] mzValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] intensityValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] m_errorValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                peakAnnotation.clear();
                for (int k = 0; k < fragmentationTable.getModel().getRowCount(); k++) {
                    mzValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 0));

                    intensityValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 1));
                    m_errorValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 2));

                    String type = (String) fragmentationTable.getModel().getValueAt(k, 3);
                    type = type.replaceFirst("frag:", "");
                    type = type.replaceFirst("ion", "");
                    type = type.replaceFirst("internal", "");

                    peakAnnotation.add(
                            new DefaultSpectrumAnnotation(
                                    mzValuesAsDouble[k],
                                    m_errorValuesAsDouble[k],
                                    Color.blue,
                                    type));
                }
                while (jGraph.getComponents().length > 0) {
                    jGraph.remove(0);
                }
                jGraph.validate();
                jGraph.repaint();
                if (jmzreader != null) {
                    try {
                        int row1 = spectrumIdentificationResultTable.getSelectedRow();
                        String sir_id = (String) spectrumIdentificationResultTable.getModel().getValueAt(row1, 0);
//                        System.out.println(sir_id);
                        SpectrumIdentificationResult spectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                        Spectrum spectrum = null;
                        String spectrumID = spectrumIdentificationResult.getSpectrumID();
                        if (sourceFile.equals("mgf")) {
                            String spectrumIndex = spectrumID.substring(6);
                            Integer index1 = Integer.valueOf(spectrumIndex) + 1;
                            spectrum = jmzreader.getSpectrumById(index1.toString());
                        }
                        if (sourceFile.equals("mzML")) {
                            spectrum = jmzreader.getSpectrumById(spectrumID);
                        }
                        peakList = spectrum.getPeakList();

                        List<Double> mzValues;
                        if (spectrum.getPeakList() != null) {
                            mzValues = new ArrayList<Double>(spectrum.getPeakList().keySet());
                        } else {
                            mzValues = Collections.emptyList();
                        }

                        double[] mz = new double[mzValues.size()];
                        double[] intensities = new double[mzValues.size()];

                        int index = 0;
                        for (Double mzValue : mzValues) {
                            mz[index] = mzValue;
                            intensities[index] = spectrum.getPeakList().get(mzValue);
                            index++;
                        }
                        spectrumPanel = new SpectrumPanel(
                                mz,
                                intensities,
                                spectrumIdentificationItem.getExperimentalMassToCharge(),
                                String.valueOf(spectrumIdentificationItem.getChargeState()),
                                spectrumIdentificationItem.getName());
                        spectrumPanel.setAnnotations(peakAnnotation);
                        jGraph.setLayout(new java.awt.BorderLayout());
                        jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                        jGraph.add(spectrumPanel);
                        jGraph.validate();
                        jGraph.repaint();
                    } catch (JMzReaderException ex) {
                        Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (JAXBException ex) {
                        Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                    }

                } else if (mzValuesAsDouble.length > 0) {
                    spectrumPanel = new SpectrumPanel(
                            mzValuesAsDouble,
                            intensityValuesAsDouble,
                            spectrumIdentificationItem.getExperimentalMassToCharge(),
                            String.valueOf(spectrumIdentificationItem.getChargeState()),
                            spectrumIdentificationItem.getName());

                    spectrumPanel.setAnnotations(peakAnnotation);

                    jGraph.setLayout(new java.awt.BorderLayout());
                    jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph.add(spectrumPanel);
                    jGraph.validate();
                    jGraph.repaint();
                    this.repaint();
                }
            }
        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    /**
     * Update Graph
     */
    private void updateGraph1() {
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTablePeptideView.getSelectedRow();
        if (row != -1) {
//            row = spectrumIdentificationItemTablePeptideView.convertRowIndexToModel(row);
            try {
                SpectrumIdentificationItem spectrumIdentificationItem = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationItem.class, (String) spectrumIdentificationItemTablePeptideView.getValueAt(row, 0));
                while (((DefaultTableModel) fragmentationTablePeptideView.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) fragmentationTablePeptideView.getModel()).removeRow(0);
                }
                for (int k = 0; k < filterListIon1.size(); k++) {
                    String nameIon = filterListIon1.get(k);
                    boolean isSelectedIon = filterCheckBoxIon1[k].isSelected();
                    for (int z = 0; z < filterListCharge1.size(); z++) {
                        String nameCharge = filterListCharge1.get(z);
                        boolean isSelectedCharge = filterCheckBoxCharge1[z].isSelected();

                        for (int i = 0; i < ionTypeList1.size(); i++) {
                            IonType ionType = ionTypeList1.get(i);
                            CvParam cvParam = ionType.getCvParam();
                            if ((nameIon.equals(cvParam.getName()) && isSelectedIon) && (nameCharge.equals(String.valueOf(ionType.getCharge())) && isSelectedCharge)) {
                                List m_mz, m_intensity, m_error;
                                m_mz = ionType.getFragmentArray().get(0).getValues();
                                m_intensity = ionType.getFragmentArray().get(1).getValues();
                                m_error = ionType.getFragmentArray().get(2).getValues();
                                if (m_mz != null && !m_mz.isEmpty()) {
                                    for (int j = 0; j < m_mz.size(); j++) {
                                        String type = cvParam.getName();
                                        ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
                                            Double.valueOf(m_mz.get(j).toString()),
                                            Double.valueOf(m_intensity.get(j).toString()),
                                            Double.valueOf(m_error.get(j).toString()),
                                            type,
                                            Integer.valueOf(ionType.getCharge())
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
//                // Start of theortical values
//                Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem.getPeptideRef());
//                String find = peptide.getPeptideSequence();
//
//                TheoreticalFragmentation tf = new TheoreticalFragmentation(find);
//
//                List<Double> tmp = new ArrayList();
//                if (bCheckBox1.isSelected() && oneCheckBox1.isSelected()) {
//                    tmp.clear();
//                    tmp = tf.getBIons(find, "1");
//                    for (int j = 0; j < tmp.size(); j++) {
//
//                        ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
//                                    Double.valueOf(tmp.get(j).toString()),
//                                    Double.valueOf("100"),
//                                    Double.valueOf("0.0"),
//                                    "T b+",
//                                    +1
//                                });
//                    }
//
//
//                }
//                if (bCheckBox1.isSelected() && twoCheckBox1.isSelected()) {
//
//                    tmp.clear();
//                    tmp = tf.getBIons(find, "2");
//                    for (int j = 0; j < tmp.size(); j++) {
//
//                        ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
//                                    Double.valueOf(tmp.get(j).toString()),
//                                    Double.valueOf("100"),
//                                    Double.valueOf("0.0"),
//                                    "T b++",
//                                    +1
//                                });
//                    }
//                }
//                if (yCheckBox1.isSelected() && oneCheckBox1.isSelected()) {
//
//
//                    tmp.clear();
//                    tmp = tf.getYIons(find, "1");
//                    for (int j = 0; j < tmp.size(); j++) {
//
//                        ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
//                                    Double.valueOf(tmp.get(j).toString()),
//                                    Double.valueOf("100"),
//                                    Double.valueOf("0.0"),
//                                    "T y+",
//                                    +1
//                                });
//                    }
//                }
//                if (yCheckBox1.isSelected() && twoCheckBox1.isSelected()) {
//
//
//                    tmp.clear();
//                    tmp = tf.getYIons(find, "2");
//                    for (int j = 0; j < tmp.size(); j++) {
//
//                        ((DefaultTableModel) fragmentationTablePeptideView.getModel()).addRow(new Object[]{
//                                    Double.valueOf(tmp.get(j).toString()),
//                                    Double.valueOf("100"),
//                                    Double.valueOf("0.0"),
//                                    "T y++",
//                                    +1
//                                });
//                    }
//                }
//
//
//
//                // End of theortical values

                peakAnnotation1.clear();
                double[] mzValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                double[] intensityValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                double[] m_errorValuesAsDouble = new double[fragmentationTablePeptideView.getModel().getRowCount()];
                peakAnnotation1.clear();
                for (int k = 0; k < fragmentationTablePeptideView.getModel().getRowCount(); k++) {
                    mzValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 0));
                    intensityValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 1));
                    m_errorValuesAsDouble[k] = (Double) (fragmentationTablePeptideView.getModel().getValueAt(k, 2));
                    String type = (String) fragmentationTablePeptideView.getModel().getValueAt(k, 3);
                    type = type.replaceFirst("frag:", "");
                    type = type.replaceFirst("ion", "");
                    type = type.replaceFirst("internal", "");
                    peakAnnotation1.add(
                            new DefaultSpectrumAnnotation(
                                    mzValuesAsDouble[k],
                                    m_errorValuesAsDouble[k],
                                    Color.blue,
                                    type));                                    // the annotation label
                }
                if (fragmentationTablePeptideView.getModel().getRowCount() > 0) {
                    spectrumPanel1 = new SpectrumPanel(
                            mzValuesAsDouble,
                            intensityValuesAsDouble,
                            spectrumIdentificationItem.getExperimentalMassToCharge(),
                            String.valueOf(spectrumIdentificationItem.getChargeState()),
                            spectrumIdentificationItem.getName());
                    spectrumPanel1.setAnnotations(peakAnnotation1);
                    while (jGraph1.getComponents().length > 0) {
                        jGraph1.remove(0);
                    }
                    jGraph1.setLayout(new java.awt.BorderLayout());
                    jGraph1.setLayout(new javax.swing.BoxLayout(jGraph1, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph1.add(spectrumPanel1);
                    jGraph1.validate();
                    jGraph1.repaint();
                    this.repaint();
                }
                setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * Update Graph
     */
    private void updateGraph() {
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTable.getSelectedRow();
        if (row != -1) {
//            row = spectrumIdentificationItemTable.convertRowIndexToModel(row);
            SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(row);
            while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
            }
            for (int k = 0; k < filterListIon.size(); k++) {
                String nameIon = filterListIon.get(k);
                boolean isSelectedIon = filterCheckBoxIon[k].isSelected();
                for (int z = 0; z < filterListCharge.size(); z++) {
                    String nameCharge = filterListCharge.get(z);
                    boolean isSelectedCharge = filterCheckBoxCharge[z].isSelected();

                    for (int i = 0; i < ionTypeList.size(); i++) {
                        IonType ionType = ionTypeList.get(i);
                        CvParam cvParam = ionType.getCvParam();
                        if ((nameIon.equals(cvParam.getName()) && isSelectedIon) && (nameCharge.equals(String.valueOf(ionType.getCharge())) && isSelectedCharge)) {
                            List m_mz, m_intensity, m_error;
                            m_mz = ionType.getFragmentArray().get(0).getValues();
                            m_intensity = ionType.getFragmentArray().get(1).getValues();
                            m_error = ionType.getFragmentArray().get(2).getValues();
                            if (m_mz != null && !m_mz.isEmpty()) {
                                for (int j = 0; j < m_mz.size(); j++) {
                                    String type = cvParam.getName();
                                    ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                        Double.valueOf(m_mz.get(j).toString()),
                                        Double.valueOf(m_intensity.get(j).toString()),
                                        Double.valueOf(m_error.get(j).toString()),
                                        type,
                                        Integer.valueOf(ionType.getCharge())
                                    });
                                }
                            }
                        }
                    }
                }
            }

            peakAnnotation.clear();
            double[] mzValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
            double[] intensityValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
            double[] m_errorValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
            peakAnnotation.clear();
            for (int k = 0; k < fragmentationTable.getModel().getRowCount(); k++) {
                mzValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 0));
                intensityValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 1));
                m_errorValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 2));
                String type = (String) fragmentationTable.getModel().getValueAt(k, 3);
                type = type.replaceFirst("frag:", "");
                type = type.replaceFirst("ion", "");
                type = type.replaceFirst("internal", "");
                peakAnnotation.add(
                        new DefaultSpectrumAnnotation(
                                mzValuesAsDouble[k],
                                m_errorValuesAsDouble[k],
                                Color.blue,
                                type));                                    // the annotation label
            }
//            if (fragmentationTable.getModel().getRowCount() > 0) {
//                spectrumPanel = new SpectrumPanel(
//                        mzValuesAsDouble,
//                        intensityValuesAsDouble,
//                        spectrumIdentificationItem.getExperimentalMassToCharge(),
//                        String.valueOf(spectrumIdentificationItem.getChargeState()),
//                        spectrumIdentificationItem.getName());
//                spectrumPanel.setAnnotations(peakAnnotation);
//                while (jGraph.getComponents().length > 0) {
//                    jGraph.remove(0);
//                }
//                jGraph.setLayout(new java.awt.BorderLayout());
//                jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
//                jGraph.add(spectrumPanel);
//                jGraph.validate();
//                jGraph.repaint();
//                this.repaint();
//            }
            while (jGraph.getComponents().length > 0) {
                jGraph.remove(0);
            }
            jGraph.validate();
            jGraph.repaint();
            if (jmzreader != null) {
                try {
                    String sir_id = (String) spectrumIdentificationResultTable.getModel().getValueAt(row, 0);
                    SpectrumIdentificationResult spectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                    Spectrum spectrum = null;
                    String spectrumID = spectrumIdentificationResult.getSpectrumID();
                    if (sourceFile.equals("mgf")) {
                        String spectrumIndex = spectrumID.substring(6);
                        Integer index1 = Integer.valueOf(spectrumIndex) + 1;
                        spectrum = jmzreader.getSpectrumById(index1.toString());
                    }
                    if (sourceFile.equals("mzML")) {
                        spectrum = jmzreader.getSpectrumById(spectrumID);
                    }
                    peakList = spectrum.getPeakList();

                    List<Double> mzValues;
                    if (spectrum.getPeakList() != null) {
                        mzValues = new ArrayList<Double>(spectrum.getPeakList().keySet());
                    } else {
                        mzValues = Collections.emptyList();
                    }

                    double[] mz = new double[mzValues.size()];
                    double[] intensities = new double[mzValues.size()];

                    int index = 0;
                    for (Double mzValue : mzValues) {
                        mz[index] = mzValue;
                        intensities[index] = spectrum.getPeakList().get(mzValue);
                        index++;
                    }
                    spectrumPanel = new SpectrumPanel(
                            mz,
                            intensities,
                            spectrumIdentificationItem.getExperimentalMassToCharge(),
                            String.valueOf(spectrumIdentificationItem.getChargeState()),
                            spectrumIdentificationItem.getName());
                    spectrumPanel.setAnnotations(peakAnnotation);
                    jGraph.setLayout(new java.awt.BorderLayout());
                    jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph.add(spectrumPanel);
                    jGraph.validate();
                    jGraph.repaint();
                } catch (JMzReaderException ex) {
                    Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                } catch (JAXBException ex) {
                    Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
                }

            } else if (mzValuesAsDouble.length > 0) {
                spectrumPanel = new SpectrumPanel(
                        mzValuesAsDouble,
                        intensityValuesAsDouble,
                        spectrumIdentificationItem.getExperimentalMassToCharge(),
                        String.valueOf(spectrumIdentificationItem.getChargeState()),
                        spectrumIdentificationItem.getName());

                spectrumPanel.setAnnotations(peakAnnotation);

                jGraph.setLayout(new java.awt.BorderLayout());
                jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                jGraph.add(spectrumPanel);
                jGraph.validate();
                jGraph.repaint();
                this.repaint();
            }
            setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

        }
    }

    /**
     * Creates spectrum Identification Item Protein Table Mouse Clicked
     */
    private void spectrumIdentificationItemProteinTableeMouseClicked(MouseEvent evt) {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

        String sii_ref = (String) spectrumIdentificationItemProteinViewTable.getValueAt(spectrumIdentificationItemProteinViewTable.getSelectedRow(), 1);

        for (int i = 0; i < spectrumIdentificationResultTable.getRowCount(); i++) {
            try {
                String sir_id = (String) spectrumIdentificationResultTable.getValueAt(i, 0);

                SpectrumIdentificationResult sir = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                List<SpectrumIdentificationItem> siiList = sir.getSpectrumIdentificationItem();
                for (int j = 0; j < siiList.size(); j++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = siiList.get(j);
                    if (sii_ref.equals(spectrumIdentificationItem.getId())) {

                        spectrumIdentificationResultTable.setRowSelectionInterval(i, i);
                        spectrumIdentificationResultTableMouseClicked(evt);

                        spectrumIdentificationItemTable.setRowSelectionInterval(j, j);
                        spectrumIdentificationItemTableMouseClicked(evt);
                        break;

                    }

                }
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }

        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        if (!secondTab) {
            loadSpectrumIdentificationResultTable();
            secondTab = true;
        }
        mainTabbedPane.setSelectedIndex(1);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fileChooser = new javax.swing.JFileChooser();
        mainPanel = new javax.swing.JPanel();
        mainTabbedPane = new javax.swing.JTabbedPane();
        proteinViewPanel = new javax.swing.JPanel();
        jProteinAmbiguityGroupPanel = new javax.swing.JPanel();
        jProteinDetectionHypothesisPanel = new javax.swing.JPanel();
        jProteinInfoPanel = new javax.swing.JPanel();
        jProteinSequencePanel = new javax.swing.JPanel();
        jProteinSequenceScrollPane = new javax.swing.JScrollPane();
        jProteinSequenceTextPane = new javax.swing.JTextPane();
        jSpectrumIdentificationItemProteinPanel = new javax.swing.JPanel();
        jScientificNameLabel = new javax.swing.JLabel();
        jScientificNameValueLabel = new javax.swing.JLabel();
        jProteinDescriptionPanel = new javax.swing.JPanel();
        jProteinDescriptionScrollPane = new javax.swing.JScrollPane();
        jProteinDescriptionEditorPane = new javax.swing.JEditorPane();
        spectrumViewPanel = new javax.swing.JPanel();
        jSplitPane1 = new javax.swing.JSplitPane();
        jSpectrumPanel = new javax.swing.JPanel();
        jFragmentationPanel = new javax.swing.JPanel();
        jGraph = new javax.swing.JPanel();
        jExperimentalFilterPanel = new javax.swing.JPanel();
        jPanel1 = new javax.swing.JPanel();
        jSpectrumIdentificationResultPanel = new javax.swing.JPanel();
        jSpectrumIdentificationItemPanel = new javax.swing.JPanel();
        jPeptideEvidencePanel = new javax.swing.JPanel();
        peptideViewPanel = new javax.swing.JPanel();
        jSplitPane2 = new javax.swing.JSplitPane();
        jSpectrumPanel1 = new javax.swing.JPanel();
        jFragmentationPanel1 = new javax.swing.JPanel();
        jGraph1 = new javax.swing.JPanel();
        jExperimentalFilterPanel1 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        psmRankLabel = new javax.swing.JLabel();
        psmRankValue = new javax.swing.JComboBox();
        jSpectrumIdentificationItemPanel1 = new javax.swing.JPanel();
        jPeptideEvidencePanel1 = new javax.swing.JPanel();
        proteinDBViewPanel = new javax.swing.JPanel();
        dBSequencePanel = new javax.swing.JPanel();
        globalStatisticsPanel = new javax.swing.JPanel();
        summaryPanel = new javax.swing.JPanel();
        totalSIRLabel = new javax.swing.JLabel();
        totalSIRLabelValue = new javax.swing.JLabel();
        totalSIILabel = new javax.swing.JLabel();
        totalSIIbelowThresholdLabel = new javax.swing.JLabel();
        totalSIIbelowThresholdLabelValue = new javax.swing.JLabel();
        totalSIIaboveThresholdLabel = new javax.swing.JLabel();
        totalSIIaboveThresholdLabelValue = new javax.swing.JLabel();
        totalSIIaboveThresholdRankOneLabel = new javax.swing.JLabel();
        totalSIIaboveThresholdRankOneLabelValue = new javax.swing.JLabel();
        percentIdentifiedSpectraLabel = new javax.swing.JLabel();
        percentIdentifiedSpectraLabelValue = new javax.swing.JLabel();
        totalPeptidesaboveThresholdLabel = new javax.swing.JLabel();
        totalPeptidesaboveThresholdLabelValue = new javax.swing.JLabel();
        totalPAGsLabel = new javax.swing.JLabel();
        totalPAGsLabelValue = new javax.swing.JLabel();
        totalPDHsLabel = new javax.swing.JLabel();
        totalPDHsaboveThresholdLabel = new javax.swing.JLabel();
        totalPDHsaboveThresholdLabelValue = new javax.swing.JLabel();
        isDecoySii = new javax.swing.JLabel();
        isDecoySiiValue = new javax.swing.JLabel();
        isDecoySiiFalse = new javax.swing.JLabel();
        isDecoySiiFalseValue = new javax.swing.JLabel();
        manualDecoyLabel = new javax.swing.JLabel();
        manualDecoyPrefix = new javax.swing.JLabel();
        manualDecoyRatio = new javax.swing.JLabel();
        manualDecoyPrefixValue = new javax.swing.JTextField();
        manualDecoyRatioValue = new javax.swing.JTextField();
        manualDecoy = new javax.swing.JCheckBox();
        fpSiiLabel = new javax.swing.JLabel();
        fpSiiValue = new javax.swing.JLabel();
        tpSiiLabel = new javax.swing.JLabel();
        tpSiiValue = new javax.swing.JLabel();
        fdrSiiLabel = new javax.swing.JLabel();
        fdrSiiValue = new javax.swing.JLabel();
        fpProteinsLabel = new javax.swing.JLabel();
        tpProteinsLabel = new javax.swing.JLabel();
        fdrProteinsLabel = new javax.swing.JLabel();
        fpProteinsValue = new javax.swing.JLabel();
        tpProteinsValue = new javax.swing.JLabel();
        fdrProteinsValue = new javax.swing.JLabel();
        manualCalculate = new javax.swing.JButton();
        fdrPanel = new javax.swing.JPanel();
        tpEvaluePanel = new javax.swing.JPanel();
        tpQvaluePanel = new javax.swing.JPanel();
        jSeparator4 = new javax.swing.JSeparator();
        siiComboBox = new javax.swing.JComboBox();
        siiListLabel = new javax.swing.JLabel();
        jComboBox1 = new javax.swing.JComboBox();
        jComboBox2 = new javax.swing.JComboBox();
        protocolPanel = new javax.swing.JPanel();
        protocolSummaryPanel = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        protocalTextPane = new javax.swing.JTextPane();
        jMenuBar = new javax.swing.JMenuBar();
        fileMenu = new javax.swing.JMenu();
        openMenuItem = new javax.swing.JMenuItem();
        saveMenuItem = new javax.swing.JMenuItem();
        jSeparator3 = new javax.swing.JPopupMenu.Separator();
        exitMenuItem = new javax.swing.JMenuItem();
        exportMenu = new javax.swing.JMenu();
        jMenuItem1 = new javax.swing.JMenuItem();
        jMenuItem2 = new javax.swing.JMenuItem();
        jMenuItem3 = new javax.swing.JMenuItem();
        exportSeparator2 = new javax.swing.JPopupMenu.Separator();
        exportFDR = new javax.swing.JMenuItem();
        optionsMenu = new javax.swing.JMenu();
        lookAndFeelMenu = new javax.swing.JMenu();
        lookAndFeelDefaultMenuItem = new javax.swing.JRadioButtonMenuItem();
        lookAndFeelMetalMenuItem = new javax.swing.JRadioButtonMenuItem();
        lookAndFeeMotifMenuItem = new javax.swing.JRadioButtonMenuItem();
        fontMenu = new javax.swing.JMenu();
        font12 = new javax.swing.JRadioButtonMenuItem();
        font14 = new javax.swing.JRadioButtonMenuItem();
        font16 = new javax.swing.JMenuItem();
        helpMenu = new javax.swing.JMenu();
        aboutMenuItem = new javax.swing.JMenuItem();

        fileChooser.setCurrentDirectory(null);
        fileChooser.addChoosableFileFilter (mzIdentMLFilter);

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("ProteoIDViewer");
        setMinimumSize(new java.awt.Dimension(900, 800));

        mainTabbedPane.setTabLayoutPolicy(javax.swing.JTabbedPane.SCROLL_TAB_LAYOUT);
        mainTabbedPane.setToolTipText("Global Statistics");
        mainTabbedPane.setPreferredSize(new java.awt.Dimension(889, 939));
        mainTabbedPane.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                mainTabbedPaneMouseClicked(evt);
            }
        });

        proteinViewPanel.setToolTipText("Protein View");
        proteinViewPanel.setName("Protein View"); // NOI18N
        proteinViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        jProteinAmbiguityGroupPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Group"));
        jProteinAmbiguityGroupPanel.setToolTipText("groups of proteins sharing some or all of the same peptides");
        jProteinAmbiguityGroupPanel.setPreferredSize(new java.awt.Dimension(772, 150));
        jProteinAmbiguityGroupPanel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jProteinAmbiguityGroupPanelMouseClicked(evt);
            }
        });

        javax.swing.GroupLayout jProteinAmbiguityGroupPanelLayout = new javax.swing.GroupLayout(jProteinAmbiguityGroupPanel);
        jProteinAmbiguityGroupPanel.setLayout(jProteinAmbiguityGroupPanelLayout);
        jProteinAmbiguityGroupPanelLayout.setHorizontalGroup(
            jProteinAmbiguityGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 614, Short.MAX_VALUE)
        );
        jProteinAmbiguityGroupPanelLayout.setVerticalGroup(
            jProteinAmbiguityGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jProteinDetectionHypothesisPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein"));
        jProteinDetectionHypothesisPanel.setToolTipText("proteins inferred based on a set of peptide spectrum matches");
        jProteinDetectionHypothesisPanel.setPreferredSize(new java.awt.Dimension(772, 150));

        javax.swing.GroupLayout jProteinDetectionHypothesisPanelLayout = new javax.swing.GroupLayout(jProteinDetectionHypothesisPanel);
        jProteinDetectionHypothesisPanel.setLayout(jProteinDetectionHypothesisPanelLayout);
        jProteinDetectionHypothesisPanelLayout.setHorizontalGroup(
            jProteinDetectionHypothesisPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 614, Short.MAX_VALUE)
        );
        jProteinDetectionHypothesisPanelLayout.setVerticalGroup(
            jProteinDetectionHypothesisPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jProteinInfoPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Info"));
        jProteinInfoPanel.setInheritsPopupMenu(true);

        jProteinSequencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Sequence"));

        jProteinSequenceTextPane.setContentType("text/html");
        jProteinSequenceTextPane.setText("");
        jProteinSequenceScrollPane.setViewportView(jProteinSequenceTextPane);

        javax.swing.GroupLayout jProteinSequencePanelLayout = new javax.swing.GroupLayout(jProteinSequencePanel);
        jProteinSequencePanel.setLayout(jProteinSequencePanelLayout);
        jProteinSequencePanelLayout.setHorizontalGroup(
            jProteinSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinSequenceScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 740, Short.MAX_VALUE)
        );
        jProteinSequencePanelLayout.setVerticalGroup(
            jProteinSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinSequenceScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 159, Short.MAX_VALUE)
        );

        jSpectrumIdentificationItemProteinPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide-Spectrum matches"));
        jSpectrumIdentificationItemProteinPanel.setToolTipText("Protein Detection Hypothesis");
        jSpectrumIdentificationItemProteinPanel.setPreferredSize(new java.awt.Dimension(772, 150));

        javax.swing.GroupLayout jSpectrumIdentificationItemProteinPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationItemProteinPanel);
        jSpectrumIdentificationItemProteinPanel.setLayout(jSpectrumIdentificationItemProteinPanelLayout);
        jSpectrumIdentificationItemProteinPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationItemProteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 740, Short.MAX_VALUE)
        );
        jSpectrumIdentificationItemProteinPanelLayout.setVerticalGroup(
            jSpectrumIdentificationItemProteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1134, Short.MAX_VALUE)
        );

        jScientificNameLabel.setText("Scientific name:");

        jProteinDescriptionPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Description"));

        jProteinDescriptionEditorPane.setContentType("text/html");
        jProteinDescriptionScrollPane.setViewportView(jProteinDescriptionEditorPane);

        javax.swing.GroupLayout jProteinDescriptionPanelLayout = new javax.swing.GroupLayout(jProteinDescriptionPanel);
        jProteinDescriptionPanel.setLayout(jProteinDescriptionPanelLayout);
        jProteinDescriptionPanelLayout.setHorizontalGroup(
            jProteinDescriptionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinDescriptionScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 740, Short.MAX_VALUE)
        );
        jProteinDescriptionPanelLayout.setVerticalGroup(
            jProteinDescriptionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinDescriptionScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 59, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jProteinInfoPanelLayout = new javax.swing.GroupLayout(jProteinInfoPanel);
        jProteinInfoPanel.setLayout(jProteinInfoPanelLayout);
        jProteinInfoPanelLayout.setHorizontalGroup(
            jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jProteinInfoPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScientificNameLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScientificNameValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 157, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(506, Short.MAX_VALUE))
            .addComponent(jProteinDescriptionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jProteinSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jSpectrumIdentificationItemProteinPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 758, Short.MAX_VALUE)
        );
        jProteinInfoPanelLayout.setVerticalGroup(
            jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jProteinInfoPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jScientificNameLabel)
                    .addComponent(jScientificNameValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 14, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(11, 11, 11)
                .addComponent(jProteinDescriptionPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jProteinSequencePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jSpectrumIdentificationItemProteinPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 1157, Short.MAX_VALUE))
        );

        jSpectrumIdentificationItemProteinPanel.getAccessibleContext().setAccessibleDescription("Spectrum Identification Item");

        javax.swing.GroupLayout proteinViewPanelLayout = new javax.swing.GroupLayout(proteinViewPanel);
        proteinViewPanel.setLayout(proteinViewPanelLayout);
        proteinViewPanelLayout.setHorizontalGroup(
            proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jProteinDetectionHypothesisPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 626, Short.MAX_VALUE)
                    .addComponent(jProteinAmbiguityGroupPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 626, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jProteinInfoPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        proteinViewPanelLayout.setVerticalGroup(
            proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, proteinViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jProteinInfoPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(proteinViewPanelLayout.createSequentialGroup()
                        .addComponent(jProteinAmbiguityGroupPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 802, Short.MAX_VALUE)
                        .addGap(18, 18, 18)
                        .addComponent(jProteinDetectionHypothesisPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 689, Short.MAX_VALUE)))
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protein View", null, proteinViewPanel, "Protein View");
        proteinViewPanel.getAccessibleContext().setAccessibleName("Protein View");
        proteinViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        spectrumViewPanel.setToolTipText("Spectrum View");
        spectrumViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        jSplitPane1.setBorder(null);
        jSplitPane1.setDividerLocation(500);

        jSpectrumPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectrum"));
        jSpectrumPanel.setAutoscrolls(true);

        jFragmentationPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Fragmentation"));
        jFragmentationPanel.setAutoscrolls(true);
        jFragmentationPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        jFragmentationPanel.setPreferredSize(new java.awt.Dimension(383, 447));

        javax.swing.GroupLayout jFragmentationPanelLayout = new javax.swing.GroupLayout(jFragmentationPanel);
        jFragmentationPanel.setLayout(jFragmentationPanelLayout);
        jFragmentationPanelLayout.setHorizontalGroup(
            jFragmentationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 893, Short.MAX_VALUE)
        );
        jFragmentationPanelLayout.setVerticalGroup(
            jFragmentationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 213, Short.MAX_VALUE)
        );

        jGraph.setBorder(javax.swing.BorderFactory.createTitledBorder("Graph"));

        javax.swing.GroupLayout jGraphLayout = new javax.swing.GroupLayout(jGraph);
        jGraph.setLayout(jGraphLayout);
        jGraphLayout.setHorizontalGroup(
            jGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jGraphLayout.setVerticalGroup(
            jGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1156, Short.MAX_VALUE)
        );

        jExperimentalFilterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Experimental Filtering"));

        javax.swing.GroupLayout jExperimentalFilterPanelLayout = new javax.swing.GroupLayout(jExperimentalFilterPanel);
        jExperimentalFilterPanel.setLayout(jExperimentalFilterPanelLayout);
        jExperimentalFilterPanelLayout.setHorizontalGroup(
            jExperimentalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jExperimentalFilterPanelLayout.setVerticalGroup(
            jExperimentalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 17, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jSpectrumPanelLayout = new javax.swing.GroupLayout(jSpectrumPanel);
        jSpectrumPanel.setLayout(jSpectrumPanelLayout);
        jSpectrumPanelLayout.setHorizontalGroup(
            jSpectrumPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jFragmentationPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 905, Short.MAX_VALUE)
            .addComponent(jExperimentalFilterPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        jSpectrumPanelLayout.setVerticalGroup(
            jSpectrumPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jSpectrumPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(24, 24, 24)
                .addComponent(jExperimentalFilterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jFragmentationPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 236, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        jSplitPane1.setRightComponent(jSpectrumPanel);

        jSpectrumIdentificationResultPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectrum List"));
        jSpectrumIdentificationResultPanel.setToolTipText("Protein Ambiguity Group");
        jSpectrumIdentificationResultPanel.setMinimumSize(new java.awt.Dimension(404, 569));

        javax.swing.GroupLayout jSpectrumIdentificationResultPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationResultPanel);
        jSpectrumIdentificationResultPanel.setLayout(jSpectrumIdentificationResultPanelLayout);
        jSpectrumIdentificationResultPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationResultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jSpectrumIdentificationResultPanelLayout.setVerticalGroup(
            jSpectrumIdentificationResultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jSpectrumIdentificationItemPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide-Spectrum matches"));
        jSpectrumIdentificationItemPanel.setToolTipText("Spectrum Identification Item");
        jSpectrumIdentificationItemPanel.setAutoscrolls(true);

        javax.swing.GroupLayout jSpectrumIdentificationItemPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationItemPanel);
        jSpectrumIdentificationItemPanel.setLayout(jSpectrumIdentificationItemPanelLayout);
        jSpectrumIdentificationItemPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationItemPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jSpectrumIdentificationItemPanelLayout.setVerticalGroup(
            jSpectrumIdentificationItemPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 224, Short.MAX_VALUE)
        );

        jPeptideEvidencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide Evidence"));
        jPeptideEvidencePanel.setToolTipText("Peptide Evidence");
        jPeptideEvidencePanel.setAutoscrolls(true);

        javax.swing.GroupLayout jPeptideEvidencePanelLayout = new javax.swing.GroupLayout(jPeptideEvidencePanel);
        jPeptideEvidencePanel.setLayout(jPeptideEvidencePanelLayout);
        jPeptideEvidencePanelLayout.setHorizontalGroup(
            jPeptideEvidencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 468, Short.MAX_VALUE)
        );
        jPeptideEvidencePanelLayout.setVerticalGroup(
            jPeptideEvidencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1029, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jPeptideEvidencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jSpectrumIdentificationResultPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jSpectrumIdentificationItemPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jSpectrumIdentificationResultPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 198, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSpectrumIdentificationItemPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPeptideEvidencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        jSpectrumIdentificationResultPanel.getAccessibleContext().setAccessibleDescription("Spectrum Identification Result");

        jSplitPane1.setLeftComponent(jPanel1);

        javax.swing.GroupLayout spectrumViewPanelLayout = new javax.swing.GroupLayout(spectrumViewPanel);
        spectrumViewPanel.setLayout(spectrumViewPanelLayout);
        spectrumViewPanelLayout.setHorizontalGroup(
            spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1)
        );
        spectrumViewPanelLayout.setVerticalGroup(
            spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1)
        );

        mainTabbedPane.addTab("Spectrum Summary", null, spectrumViewPanel, "Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleName("Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleDescription("Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        peptideViewPanel.setToolTipText("Peptide View");
        peptideViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        jSplitPane2.setBorder(null);
        jSplitPane2.setDividerLocation(500);

        jSpectrumPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectrum"));
        jSpectrumPanel1.setAutoscrolls(true);
        jSpectrumPanel1.setPreferredSize(new java.awt.Dimension(362, 569));

        jFragmentationPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Fragmentation"));
        jFragmentationPanel1.setAutoscrolls(true);
        jFragmentationPanel1.setPreferredSize(new java.awt.Dimension(383, 447));

        javax.swing.GroupLayout jFragmentationPanel1Layout = new javax.swing.GroupLayout(jFragmentationPanel1);
        jFragmentationPanel1.setLayout(jFragmentationPanel1Layout);
        jFragmentationPanel1Layout.setHorizontalGroup(
            jFragmentationPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 893, Short.MAX_VALUE)
        );
        jFragmentationPanel1Layout.setVerticalGroup(
            jFragmentationPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 227, Short.MAX_VALUE)
        );

        jGraph1.setBorder(javax.swing.BorderFactory.createTitledBorder("Graph"));

        javax.swing.GroupLayout jGraph1Layout = new javax.swing.GroupLayout(jGraph1);
        jGraph1.setLayout(jGraph1Layout);
        jGraph1Layout.setHorizontalGroup(
            jGraph1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jGraph1Layout.setVerticalGroup(
            jGraph1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1160, Short.MAX_VALUE)
        );

        jExperimentalFilterPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Experimental Filtering"));

        javax.swing.GroupLayout jExperimentalFilterPanel1Layout = new javax.swing.GroupLayout(jExperimentalFilterPanel1);
        jExperimentalFilterPanel1.setLayout(jExperimentalFilterPanel1Layout);
        jExperimentalFilterPanel1Layout.setHorizontalGroup(
            jExperimentalFilterPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 894, Short.MAX_VALUE)
        );
        jExperimentalFilterPanel1Layout.setVerticalGroup(
            jExperimentalFilterPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 17, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jSpectrumPanel1Layout = new javax.swing.GroupLayout(jSpectrumPanel1);
        jSpectrumPanel1.setLayout(jSpectrumPanel1Layout);
        jSpectrumPanel1Layout.setHorizontalGroup(
            jSpectrumPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jExperimentalFilterPanel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jFragmentationPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, 905, Short.MAX_VALUE)
            .addComponent(jGraph1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        jSpectrumPanel1Layout.setVerticalGroup(
            jSpectrumPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jSpectrumPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jGraph1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jExperimentalFilterPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jFragmentationPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        jSplitPane2.setRightComponent(jSpectrumPanel1);

        psmRankLabel.setText("Peptide-Spectrum matches with Rank: ");

        psmRankValue.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "<=1", "<=2", "<=3", "All" }));
        psmRankValue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                psmRankValueActionPerformed(evt);
            }
        });

        jSpectrumIdentificationItemPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide-Spectrum matches"));
        jSpectrumIdentificationItemPanel1.setToolTipText("Spectrum Identification Item");
        jSpectrumIdentificationItemPanel1.setAutoscrolls(true);

        javax.swing.GroupLayout jSpectrumIdentificationItemPanel1Layout = new javax.swing.GroupLayout(jSpectrumIdentificationItemPanel1);
        jSpectrumIdentificationItemPanel1.setLayout(jSpectrumIdentificationItemPanel1Layout);
        jSpectrumIdentificationItemPanel1Layout.setHorizontalGroup(
            jSpectrumIdentificationItemPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jSpectrumIdentificationItemPanel1Layout.setVerticalGroup(
            jSpectrumIdentificationItemPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 303, Short.MAX_VALUE)
        );

        jPeptideEvidencePanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide Evidence"));
        jPeptideEvidencePanel1.setToolTipText("Peptide Evidence");
        jPeptideEvidencePanel1.setAutoscrolls(true);

        javax.swing.GroupLayout jPeptideEvidencePanel1Layout = new javax.swing.GroupLayout(jPeptideEvidencePanel1);
        jPeptideEvidencePanel1.setLayout(jPeptideEvidencePanel1Layout);
        jPeptideEvidencePanel1Layout.setHorizontalGroup(
            jPeptideEvidencePanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jPeptideEvidencePanel1Layout.setVerticalGroup(
            jPeptideEvidencePanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1124, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addComponent(psmRankLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(psmRankValue, javax.swing.GroupLayout.PREFERRED_SIZE, 79, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 212, Short.MAX_VALUE))
                    .addComponent(jSpectrumIdentificationItemPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPeptideEvidencePanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGap(14, 14, 14)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(psmRankLabel)
                    .addComponent(psmRankValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSpectrumIdentificationItemPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jPeptideEvidencePanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jSplitPane2.setLeftComponent(jPanel2);

        javax.swing.GroupLayout peptideViewPanelLayout = new javax.swing.GroupLayout(peptideViewPanel);
        peptideViewPanel.setLayout(peptideViewPanelLayout);
        peptideViewPanelLayout.setHorizontalGroup(
            peptideViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 1422, Short.MAX_VALUE)
        );
        peptideViewPanelLayout.setVerticalGroup(
            peptideViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane2)
        );

        mainTabbedPane.addTab("Peptide Summary", null, peptideViewPanel, "Peptide Summary");

        proteinDBViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        dBSequencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("DB Sequence"));

        javax.swing.GroupLayout dBSequencePanelLayout = new javax.swing.GroupLayout(dBSequencePanel);
        dBSequencePanel.setLayout(dBSequencePanelLayout);
        dBSequencePanelLayout.setHorizontalGroup(
            dBSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1390, Short.MAX_VALUE)
        );
        dBSequencePanelLayout.setVerticalGroup(
            dBSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1486, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout proteinDBViewPanelLayout = new javax.swing.GroupLayout(proteinDBViewPanel);
        proteinDBViewPanel.setLayout(proteinDBViewPanelLayout);
        proteinDBViewPanelLayout.setHorizontalGroup(
            proteinDBViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinDBViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(dBSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        proteinDBViewPanelLayout.setVerticalGroup(
            proteinDBViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinDBViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(dBSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protein DB View", proteinDBViewPanel);
        proteinDBViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        globalStatisticsPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        summaryPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Summary"));

        totalSIRLabel.setText("Total SIR:");

        totalSIRLabelValue.setText("0");

        totalSIILabel.setText("Total PSM:");

        totalSIIbelowThresholdLabel.setText("Total PSM below Threshold:");

        totalSIIbelowThresholdLabelValue.setText("0");

        totalSIIaboveThresholdLabel.setText("Total PSM pass Threshold:");

        totalSIIaboveThresholdLabelValue.setText("0");

        totalSIIaboveThresholdRankOneLabel.setText("Total PSM pass Threshold and where rank=1:");

        totalSIIaboveThresholdRankOneLabelValue.setText("0");

        percentIdentifiedSpectraLabel.setText("Percent identified spectra:");

        percentIdentifiedSpectraLabelValue.setText("0");

        totalPeptidesaboveThresholdLabel.setText("Total non-redundant peptides above Threshold:");

        totalPeptidesaboveThresholdLabelValue.setText("0");

        totalPAGsLabel.setText("Total PAGs:");

        totalPAGsLabelValue.setText("0");

        totalPDHsLabel.setText("Total PDHs:");

        totalPDHsaboveThresholdLabel.setText("Total PDH above Threshold:");

        totalPDHsaboveThresholdLabelValue.setText("0");

        isDecoySii.setText("PSM with Decoy = true:");

        isDecoySiiValue.setText("0");

        isDecoySiiFalse.setText("PSM with Decoy = false:");

        isDecoySiiFalseValue.setText("0");

        manualDecoyLabel.setText("Manual Decoy:");

        manualDecoyPrefix.setText("Prefix");

        manualDecoyRatio.setText("Ratio");

        manualDecoyPrefixValue.setText("Rev_");
        manualDecoyPrefixValue.setEnabled(false);

        manualDecoyRatioValue.setText("1");
        manualDecoyRatioValue.setEnabled(false);

        manualDecoy.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                manualDecoyActionPerformed(evt);
            }
        });

        fpSiiLabel.setText("FP for PSM:");

        fpSiiValue.setText("0");

        tpSiiLabel.setText("TP for PSM:");

        tpSiiValue.setText("0");

        fdrSiiLabel.setText("FDR for PSM:");

        fdrSiiValue.setText("0");

        fpProteinsLabel.setText("FP for proteins:");

        tpProteinsLabel.setText("TP for proteins:");

        fdrProteinsLabel.setText("FDR for proteins:");

        fpProteinsValue.setText("0");

        tpProteinsValue.setText("0");

        fdrProteinsValue.setText("0");

        manualCalculate.setText("Calculate / Show graphs");
        manualCalculate.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                manualCalculateActionPerformed(evt);
            }
        });

        fdrPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("FDR Graph"));
        fdrPanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout fdrPanelLayout = new javax.swing.GroupLayout(fdrPanel);
        fdrPanel.setLayout(fdrPanelLayout);
        fdrPanelLayout.setHorizontalGroup(
            fdrPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        fdrPanelLayout.setVerticalGroup(
            fdrPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        tpEvaluePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("TP vs FP vs E-value"));
        tpEvaluePanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout tpEvaluePanelLayout = new javax.swing.GroupLayout(tpEvaluePanel);
        tpEvaluePanel.setLayout(tpEvaluePanelLayout);
        tpEvaluePanelLayout.setHorizontalGroup(
            tpEvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        tpEvaluePanelLayout.setVerticalGroup(
            tpEvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 917, Short.MAX_VALUE)
        );

        tpQvaluePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("TP vs Q-value"));
        tpQvaluePanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout tpQvaluePanelLayout = new javax.swing.GroupLayout(tpQvaluePanel);
        tpQvaluePanel.setLayout(tpQvaluePanelLayout);
        tpQvaluePanelLayout.setHorizontalGroup(
            tpQvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        tpQvaluePanelLayout.setVerticalGroup(
            tpQvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 918, Short.MAX_VALUE)
        );

        jSeparator4.setOrientation(javax.swing.SwingConstants.VERTICAL);

        siiComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                siiComboBoxActionPerformed(evt);
            }
        });

        siiListLabel.setText("SII List:");

        jComboBox2.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Better scores are lower", "Better scores are higher" }));

        javax.swing.GroupLayout summaryPanelLayout = new javax.swing.GroupLayout(summaryPanel);
        summaryPanel.setLayout(summaryPanelLayout);
        summaryPanelLayout.setHorizontalGroup(
            summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(summaryPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addComponent(fdrPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 450, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(tpEvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 450, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(tpQvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 450, Short.MAX_VALUE)
                        .addContainerGap())
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(percentIdentifiedSpectraLabel)
                            .addComponent(totalPeptidesaboveThresholdLabel)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(totalPAGsLabel)
                                    .addComponent(totalPDHsLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(totalPDHsaboveThresholdLabel)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 408, Short.MAX_VALUE)
                        .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(isDecoySii)
                                    .addComponent(isDecoySiiFalse))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(isDecoySiiFalseValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(isDecoySiiValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(fpSiiLabel)
                                    .addComponent(tpSiiLabel)
                                    .addComponent(fdrSiiLabel)
                                    .addComponent(fpProteinsLabel)
                                    .addComponent(tpProteinsLabel)
                                    .addComponent(fdrProteinsLabel))
                                .addGap(40, 40, 40)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(fpSiiValue, javax.swing.GroupLayout.PREFERRED_SIZE, 39, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(tpSiiValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(fpProteinsValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(tpProteinsValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(fdrProteinsValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(fdrSiiValue, javax.swing.GroupLayout.PREFERRED_SIZE, 47, javax.swing.GroupLayout.PREFERRED_SIZE))))
                        .addGap(552, 552, 552))
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addComponent(manualDecoy)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(manualDecoyLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(manualDecoyPrefix)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(manualDecoyPrefixValue, javax.swing.GroupLayout.PREFERRED_SIZE, 51, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(totalSIRLabel)
                                    .addComponent(totalSIILabel)
                                    .addComponent(siiListLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(siiComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(totalSIRLabelValue)
                                            .addComponent(totalSIIbelowThresholdLabel)
                                            .addComponent(totalSIIaboveThresholdLabel)
                                            .addComponent(totalSIIaboveThresholdRankOneLabel))
                                        .addGap(22, 22, 22)
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                            .addComponent(totalSIIbelowThresholdLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                                            .addComponent(totalSIIaboveThresholdLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(totalSIIaboveThresholdRankOneLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(totalPeptidesaboveThresholdLabelValue, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(percentIdentifiedSpectraLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(totalPDHsaboveThresholdLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(totalPAGsLabelValue, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addGap(150, 150, 150)
                                        .addComponent(manualDecoyRatio, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(manualDecoyRatioValue, javax.swing.GroupLayout.PREFERRED_SIZE, 51, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(jComboBox1, javax.swing.GroupLayout.PREFERRED_SIZE, 155, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jComboBox2, javax.swing.GroupLayout.PREFERRED_SIZE, 177, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(manualCalculate, javax.swing.GroupLayout.PREFERRED_SIZE, 162, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addContainerGap())))
        );
        summaryPanelLayout.setVerticalGroup(
            summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(summaryPanelLayout.createSequentialGroup()
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(siiComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(siiListLabel))
                                .addGap(8, 8, 8)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIRLabel)
                                    .addComponent(totalSIRLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIILabel)
                                    .addComponent(totalSIIbelowThresholdLabel)
                                    .addComponent(totalSIIbelowThresholdLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIIaboveThresholdLabel)
                                    .addComponent(totalSIIaboveThresholdLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIIaboveThresholdRankOneLabel)
                                    .addComponent(totalSIIaboveThresholdRankOneLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addComponent(totalPAGsLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                            .addComponent(totalPDHsLabel)
                                            .addComponent(totalPDHsaboveThresholdLabel)))
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addComponent(totalPAGsLabelValue)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(totalPDHsaboveThresholdLabelValue)))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(percentIdentifiedSpectraLabelValue)
                                    .addComponent(percentIdentifiedSpectraLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalPeptidesaboveThresholdLabel)
                                    .addComponent(totalPeptidesaboveThresholdLabelValue)))
                            .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, 237, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addGap(22, 22, 22)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(isDecoySii)
                            .addComponent(isDecoySiiValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(isDecoySiiFalse)
                            .addComponent(isDecoySiiFalseValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(fpSiiLabel)
                            .addComponent(fpSiiValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(tpSiiLabel)
                            .addComponent(tpSiiValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(fdrSiiLabel)
                            .addComponent(fdrSiiValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(fpProteinsLabel)
                            .addComponent(fpProteinsValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(tpProteinsLabel)
                            .addComponent(tpProteinsValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(fdrProteinsLabel)
                            .addComponent(fdrProteinsValue))))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 253, Short.MAX_VALUE)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(manualDecoyPrefix)
                        .addComponent(manualDecoyPrefixValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(manualDecoyLabel)
                        .addComponent(manualDecoyRatio)
                        .addComponent(manualDecoyRatioValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(manualCalculate)
                        .addComponent(jComboBox1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jComboBox2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(manualDecoy))
                .addGap(11, 11, 11)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(tpEvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 940, Short.MAX_VALUE)
                    .addComponent(tpQvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 940, Short.MAX_VALUE)
                    .addComponent(fdrPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 940, Short.MAX_VALUE))
                .addContainerGap())
        );

        javax.swing.GroupLayout globalStatisticsPanelLayout = new javax.swing.GroupLayout(globalStatisticsPanel);
        globalStatisticsPanel.setLayout(globalStatisticsPanelLayout);
        globalStatisticsPanelLayout.setHorizontalGroup(
            globalStatisticsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(globalStatisticsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(summaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        globalStatisticsPanelLayout.setVerticalGroup(
            globalStatisticsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(globalStatisticsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(summaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Global Statistics", null, globalStatisticsPanel, "Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleName("Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleDescription("Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        protocolPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        protocolSummaryPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Summary"));

        protocalTextPane.setContentType("text/html");
        protocalTextPane.setEditable(false);
        jScrollPane1.setViewportView(protocalTextPane);

        javax.swing.GroupLayout protocolSummaryPanelLayout = new javax.swing.GroupLayout(protocolSummaryPanel);
        protocolSummaryPanel.setLayout(protocolSummaryPanelLayout);
        protocolSummaryPanelLayout.setHorizontalGroup(
            protocolSummaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolSummaryPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 1370, Short.MAX_VALUE)
                .addContainerGap())
        );
        protocolSummaryPanelLayout.setVerticalGroup(
            protocolSummaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolSummaryPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 1464, Short.MAX_VALUE)
                .addContainerGap())
        );

        javax.swing.GroupLayout protocolPanelLayout = new javax.swing.GroupLayout(protocolPanel);
        protocolPanel.setLayout(protocolPanelLayout);
        protocolPanelLayout.setHorizontalGroup(
            protocolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(protocolSummaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        protocolPanelLayout.setVerticalGroup(
            protocolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(protocolSummaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protocols", protocolPanel);
        protocolPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        javax.swing.GroupLayout mainPanelLayout = new javax.swing.GroupLayout(mainPanel);
        mainPanel.setLayout(mainPanelLayout);
        mainPanelLayout.setHorizontalGroup(
            mainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainPanelLayout.createSequentialGroup()
                .addComponent(mainTabbedPane, javax.swing.GroupLayout.DEFAULT_SIZE, 1427, Short.MAX_VALUE)
                .addContainerGap())
        );
        mainPanelLayout.setVerticalGroup(
            mainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainPanelLayout.createSequentialGroup()
                .addComponent(mainTabbedPane, javax.swing.GroupLayout.DEFAULT_SIZE, 1559, Short.MAX_VALUE)
                .addContainerGap())
        );

        fileMenu.setText("File");

        openMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_O, java.awt.event.InputEvent.SHIFT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        openMenuItem.setText("Open");
        openMenuItem.setToolTipText("Open");
        openMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(openMenuItem);

        saveMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.SHIFT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        saveMenuItem.setText("Save");
        saveMenuItem.setEnabled(false);
        saveMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(saveMenuItem);
        fileMenu.add(jSeparator3);

        exitMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_Q, java.awt.event.InputEvent.CTRL_MASK));
        exitMenuItem.setText("Exit");
        exitMenuItem.setToolTipText("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(exitMenuItem);

        jMenuBar.add(fileMenu);

        exportMenu.setText("Export");

        jMenuItem1.setText("Export Proteins Only");
        jMenuItem1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem1ActionPerformed(evt);
            }
        });
        exportMenu.add(jMenuItem1);

        jMenuItem2.setText("Export Protein Groups");
        jMenuItem2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem2ActionPerformed(evt);
            }
        });
        exportMenu.add(jMenuItem2);

        jMenuItem3.setText("Export PSMs");
        jMenuItem3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem3ActionPerformed(evt);
            }
        });
        exportMenu.add(jMenuItem3);
        exportMenu.add(exportSeparator2);

        exportFDR.setText("Export FDR as CSV");
        exportFDR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportFDRActionPerformed(evt);
            }
        });
        exportMenu.add(exportFDR);

        jMenuBar.add(exportMenu);

        optionsMenu.setText("Options");

        lookAndFeelMenu.setText("Look & Feel");

        lookAndFeelDefaultMenuItem.setSelected(true);
        lookAndFeelDefaultMenuItem.setText("Default");
        lookAndFeelDefaultMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeelDefaultMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeelDefaultMenuItem);

        lookAndFeelMetalMenuItem.setText("Metal");
        lookAndFeelMetalMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeelMetalMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeelMetalMenuItem);

        lookAndFeeMotifMenuItem.setText("Motif");
        lookAndFeeMotifMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeeMotifMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeeMotifMenuItem);

        optionsMenu.add(lookAndFeelMenu);

        fontMenu.setText("Font");

        font12.setText("12");
        font12.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font12ActionPerformed(evt);
            }
        });
        fontMenu.add(font12);

        font14.setText("14");
        font14.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font14ActionPerformed(evt);
            }
        });
        fontMenu.add(font14);

        font16.setText("16");
        font16.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font16ActionPerformed(evt);
            }
        });
        fontMenu.add(font16);

        optionsMenu.add(fontMenu);

        jMenuBar.add(optionsMenu);

        helpMenu.setText("Help");

        aboutMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_A, java.awt.event.InputEvent.CTRL_MASK));
        aboutMenuItem.setText("About");
        aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                aboutMenuItemActionPerformed(evt);
            }
        });
        helpMenu.add(aboutMenuItem);

        jMenuBar.add(helpMenu);

        setJMenuBar(jMenuBar);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(mainPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(mainPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        getAccessibleContext().setAccessibleDescription("MzIdentML Viewer");

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * exit Menu Item Action Performed
     */
    private void exitMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitMenuItemActionPerformed
        int exit = JOptionPane.showConfirmDialog(this, "Are you sure?", "Close ProteoIDViewer", JOptionPane.YES_NO_OPTION);
        if (exit == JOptionPane.YES_OPTION) {
            System.exit(0);
        }
    }//GEN-LAST:event_exitMenuItemActionPerformed
    /**
     * about Menu Item Action Performed
     */
    private void aboutMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_aboutMenuItemActionPerformed
        JOptionPane.showMessageDialog(this, "ProteoIDViewer \n\nPost-genomic bioinformatics group\nInstitute of Integrative Biology\nUniversity of Liverpool\n", "About", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_aboutMenuItemActionPerformed
    /**
     * protein Detection Hypothesis Table Mouse Clicked
     */
    private void proteinDetectionHypothesisTableMouseClicked(MouseEvent evt) {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = proteinDetectionHypothesisTable.getSelectedRow();
        SpectrumIdentificationItem spectrumIdentificationItem2 = null;
        if (row != -1) {
//            row = spectrumIdentificationItemTable.convertRowIndexToModel(row);
            try {
                while (spectrumIdentificationItemProteinViewTable.getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
                }
                spectrumIdentificationItemProteinViewTable.scrollRowToVisible(0);
//                row = proteinDetectionHypothesisTable.convertRowIndexToModel(row);
                ProteinDetectionHypothesis proteinDetectionHypothesis = mzIdentMLUnmarshaller.unmarshal(ProteinDetectionHypothesis.class, (String) proteinDetectionHypothesisTable.getModel().getValueAt(row, 0));
//                System.out.println((String) proteinDetectionHypothesisTable.getModel().getValueAt(row, 0));
                DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, proteinDetectionHypothesis.getDBSequenceRef());
//                System.out.println(proteinDetectionHypothesis.getDBSequenceRef());
//                System.out.println(dBSequence.getAccession());
                List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesis.getPeptideHypothesis();
                String proteinSequence = "";
                String protein_description = "";
                if (dBSequence != null) {
                    List<CvParam> cvParamListDBSequence = dBSequence.getCvParam();
                    String scientific_name = null;

                    for (int j = 0; j < cvParamListDBSequence.size(); j++) {
                        CvParam cvParam = cvParamListDBSequence.get(j);
                        if (cvParam.getName().equals("taxonomy: scientific name")) {
                            scientific_name = cvParam.getValue();
                        }
                        if (cvParam.getName().equals("protein description")) {
                            protein_description = cvParam.getValue();
                        }
                    }
                    jScientificNameValueLabel.setText(scientific_name);
                    jProteinDescriptionEditorPane.setText(protein_description);
                    proteinSequence = dBSequence.getSeq();

                    for (int i = 0; i < peptideHypothesisList.size(); i++) {
                        PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(i);
                        for (int j = 0; j < peptideHypothesis.getSpectrumIdentificationItemRef().size(); j++) {
                            spectrumIdentificationItem2 = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationItem.class, peptideHypothesis.getSpectrumIdentificationItemRef().get(j).getSpectrumIdentificationItemRef());
                        }
                        List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem2.getPeptideEvidenceRef();

                        for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);
                            PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                            if (peptideEvidence.getDBSequenceRef().equals(proteinDetectionHypothesis.getDBSequenceRef())) {
                                Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem2.getPeptideRef());
                                if (peptide != null) {
                                    List<Modification> modificationList = peptide.getModification();
                                    Modification modification;
                                    String residues = null;
                                    Integer location = null;
                                    String modificationName = null;
                                    CvParam modificationCvParam;
                                    String combine = null;
                                    if (modificationList.size() > 0) {
                                        modification = modificationList.get(0);
                                        location = modification.getLocation();
                                        if (modification.getResidues().size() > 0) {
                                            residues = modification.getResidues().get(0);
                                        }
                                        List<CvParam> modificationCvParamList = modification.getCvParam();
                                        if (modificationCvParamList.size() > 0) {
                                            modificationCvParam = modificationCvParamList.get(0);
                                            modificationName = modificationCvParam.getName();
                                        }
                                    }
                                    if (modificationName != null) {
                                        combine = modificationName;
                                    }
                                    if (residues != null) {
                                        combine = combine + " on residues: " + residues;
                                    }
                                    if (location != null) {
                                        combine = combine + " at location: " + location;
                                    }

                                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).addRow(new String[]{
                                        peptide.getPeptideSequence(), spectrumIdentificationItem2.getId(), combine});
                                    String find = peptide.getPeptideSequence();
//                                    String replace = "<FONT COLOR=\"red\">" + find + "</FONT>";
                                    Pattern pattern = Pattern.compile(find);
                                    if (proteinSequence != null) {
                                        Matcher matcher = pattern.matcher(proteinSequence);
//                                        proteinSequence = matcher.replaceAll(replace);
                                    }
                                    List<CvParam> cvParamListSpectrumIdentificationItem = spectrumIdentificationItem2.getCvParam();

                                    for (int s = 0; s < cvParamListSpectrumIdentificationItem.size(); s++) {
                                        CvParam cvParam = cvParamListSpectrumIdentificationItem.get(s);
                                        String newCol = cvParam.getName();
                                        if (newCol.equals("mascot:score")) {
                                            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).setValueAt(Double.valueOf(cvParam.getValue()), ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() - 1, 3);
                                        }
                                        String accession = cvParam.getAccession();
                                        if (accession.equals("MS:1001330")
                                                || accession.equals("MS:1001172")
                                                || accession.equals("MS:1001159")
                                                || accession.equals("MS:1001328")) {
                                            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).setValueAt(Double.valueOf(cvParam.getValue()), ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() - 1, 4);
                                        }

                                    }
                                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).setValueAt(spectrumIdentificationItem2.isPassThreshold(), ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() - 1, 5);
                                }

                            }

                        }
                    }
                }
                if (proteinSequence != null) {
                    StringBuilder sb = new StringBuilder(proteinSequence);
                    StringBuilder sb_new = new StringBuilder();
                    int i = 0;

                    for (int j = 0; j < sb.length(); j++) {

                        if (i % 60 == 0 && i != 0) {
                            sb_new.append("<BR>");
                        }
                        i = i + 1;
                        sb_new.append(sb.charAt(j));
//                        if (sb.charAt(j) == '<') {
//                            if (sb.charAt(j + 1) == '/') {
//                                sb_new.append(sb.charAt(j + 1));
//                                for (int z = j + 2; z <= j + 7; z++) {
//                                    sb_new.append(sb.charAt(z));
//                                }
//                                j = j + 7;
//                            } else {
//                                for (int z = j + 1; z <= j + 18; z++) {
//                                    sb_new.append(sb.charAt(z));
//                                }
//                                j = j + 18;
//                            }
//                        }

//                        if (i % 60 == 0 && i != 0) {
//                            sb_new.append("<BR>");
//                        }
//                        i = i + 1;
//                        if (sb.charAt(j) == '<') {
//                         while (sb.charAt(j+1)!='>')   {
//                             j=j+1;
//                        }
//                        sb_new.append(sb.charAt(j));
//                       }
                    }
                    jProteinSequenceTextPane.setText("<FONT FACE=\"Courier New\">" + sb_new.toString() + "</FONT>");
//                    System.out.println("-----------------------");
//                    System.out.println(jProteinSequenceTextPane.getText());

                }
//                SimpleAttributeSet sa = new SimpleAttributeSet();
//                StyleConstants.setAlignment(sa, StyleConstants.ALIGN_JUSTIFIED);
//
//                jProteinSequenceTextPane.getStyledDocument().setParagraphAttributes(0, 60, sa, false);
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    /**
     * protein Ambiguity Group Table Mouse Clicked
     */
    private void proteinAmbiguityGroupTableMouseClicked(MouseEvent evt) {
        jProteinSequenceTextPane.setText("");
        jProteinSequenceTextPane.setText("");
        jScientificNameValueLabel.setText("");
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = proteinAmbiguityGroupTable.getSelectedRow();

        if (row != -1) {
            row = proteinAmbiguityGroupTable.convertRowIndexToModel(row);
            try {
                while (proteinDetectionHypothesisTable.getRowCount() > 0) {
                    ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
                }
                proteinDetectionHypothesisTable.scrollRowToVisible(0);
                String pag_id = (String) proteinAmbiguityGroupTable.getModel().getValueAt(row, 0);

                ProteinAmbiguityGroup proteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshal(ProteinAmbiguityGroup.class, pag_id);

                List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
                if (proteinDetectionHypothesisList.size() > 0) {
                    for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {

                        ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
                        DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, proteinDetectionHypothesis.getDBSequenceRef());
                        boolean isDecoy = checkIfProteinDetectionHypothesisIsDecoy(proteinDetectionHypothesis);
                        List<CvParam> cvParamList = proteinDetectionHypothesis.getCvParam();
                        String score = " ";
                        String number_peptide = " ";
                        for (int j = 0; j < cvParamList.size(); j++) {
                            CvParam cvParam = cvParamList.get(j);
                            if (cvParam.getName().contains("score")) {
                                score = cvParam.getValue();
                            }

                        }
                        String dBSequenceAccession = "";
                        if (dBSequence != null) {
                            dBSequenceAccession = dBSequence.getAccession();
                        }
                        if (proteinDetectionHypothesis.getPeptideHypothesis() != null) {
                            number_peptide = String.valueOf(proteinDetectionHypothesis.getPeptideHypothesis().size());
                        }
                        
                   
                        ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).addRow(new Object[]{
                            proteinDetectionHypothesis.getId(),
                            dBSequenceAccession,
                            score,
                            "",
                            Integer.valueOf(number_peptide),
                            isDecoy,
                            proteinDetectionHypothesis.isPassThreshold()
                        });
                    }
                }
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private boolean checkIfProteinDetectionHypothesisIsDecoy(ProteinDetectionHypothesis proteinDetectionHypothesis) {
        boolean result = false;
        List<PeptideHypothesis> PeptideHyposthesisList = proteinDetectionHypothesis.getPeptideHypothesis();
        for (int i = 0; i < PeptideHyposthesisList.size(); i++) {
            try {
                PeptideHypothesis peptideHypothesis = PeptideHyposthesisList.get(i);
                String peptidRef = peptideHypothesis.getPeptideEvidenceRef();
                PeptideEvidence peptiedEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptidRef);
                if (peptiedEvidence.isIsDecoy()) {
                    result = true;
                    break;
                }
            } catch (JAXBException ex) {
                Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        return result;
    }

    public void openCLI(String filename) {

        progressBarDialog = new ProgressBarDialog(this, true);
        final Thread thread = new Thread(new Runnable() {

            @Override
            public void run() {
                progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                progressBarDialog.setVisible(true);
            }
        }, "ProgressBarDialog");

        final File mzid_file = new File(filename);
        setTitle("ProteoIDViewer   -  " + mzid_file.getPath());
        thread.start();

        new Thread("LoadingThread") {

            @Override
            public void run() {
                try {
                    if (mzid_file.getPath().endsWith(".gz")) {
                        File outFile = null;
                        FileOutputStream fos = null;

                        GZIPInputStream gin = new GZIPInputStream(new FileInputStream(mzid_file));
                        outFile = new File(mzid_file.getParent(), mzid_file.getName().replaceAll("\\.gz$", ""));
                        fos = new FileOutputStream(outFile);
                        byte[] buf = new byte[100000];
                        int len;
                        while ((len = gin.read(buf)) > 0) {
                            fos.write(buf, 0, len);
                        }
                        fos.close();

                        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                        fileName = outFile.getAbsolutePath();
                    } else if (mzid_file.getPath().endsWith(".omx")) {
                        File outFile = null;
                        outFile = new File(fileChooser.getCurrentDirectory(), mzid_file.getName().replaceAll(".omx", ".mzid"));
                        new Omssa2mzid(mzid_file.getPath(), outFile.getPath(), false);
                        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                        fileName = outFile.getAbsolutePath();
                    } else if (mzid_file.getPath().endsWith(".xml")) {
                        File outFile = null;
                        outFile = new File(fileChooser.getCurrentDirectory(), mzid_file.getName().replaceAll(".omx", ".mzid"));
                        new Tandem2mzid(mzid_file.getPath(), outFile.getPath());
                        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                        fileName = outFile.getAbsolutePath();
                    } else {
                        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(mzid_file);
                        fileName = mzid_file.getAbsolutePath();
                    }

                    if (!mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.1.") && !mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.")) {
                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        JOptionPane.showMessageDialog(null, "The file is not compatible with the Viewer: different mzIdentMl version", "mzIdentMl version", JOptionPane.INFORMATION_MESSAGE);
                        return;
                    }
                    jmzreader = null;
                    createTables();
                    clearSummaryStats();
                    mainTabbedPane.setSelectedIndex(0);
                    secondTab = false;
                    thirdTab = false;
                    fourthTab = false;
                    fifthTab = false;
                    sixthTab = false;
                    loadProteinAmbiguityGroupTable();

                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();

                    String message = "Do you want to load spectrum source file?";

                    int answer = JOptionPane.showConfirmDialog(null, message);
                    if (answer == JOptionPane.YES_OPTION) {
                        JFileChooser fc;
                        //Create a file chooser
                        fc = new JFileChooser();
                        fc.setCurrentDirectory(fileChooser.getCurrentDirectory());

                        fc.addChoosableFileFilter(new SourceFileFilter());
                        int returnVal1 = fc.showOpenDialog(null);

                        if (returnVal1 == JFileChooser.APPROVE_OPTION) {
                            try {
                                File file = fc.getSelectedFile();
                                if (file.getAbsolutePath().toLowerCase().endsWith("mgf")) {
                                    jmzreader = new MgfFile(file);
                                    sourceFile = "mgf";
                                    JOptionPane.showMessageDialog(null, file.getName() + " is loaded", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                } else if (file.getAbsolutePath().toLowerCase().endsWith("mzml")) {
                                    jmzreader = new MzMlWrapper(file);
                                    sourceFile = "mzML";
                                    JOptionPane.showMessageDialog(null, file.getName() + " is loaded", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                } else {
                                    JOptionPane.showMessageDialog(null, file.getName() + " is not supported", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                }
                            } catch (JMzReaderException ex) {
                                System.out.println(ex.getMessage());
                            }
                        }
                    }

                    if (proteinAmbiguityGroupTable.getRowCount() == 0) {
                        JOptionPane.showMessageDialog(null, "There is no protein view for this file", "Protein View", JOptionPane.INFORMATION_MESSAGE);
                    }
                    //loadSummaryStats();

                } catch (OutOfMemoryError error) {
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                    Runtime.getRuntime().gc();
                    JOptionPane.showMessageDialog(null, "Out of Memory Error.", "Error", JOptionPane.ERROR_MESSAGE);

                    System.exit(0);
                } catch (Exception ex) {
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                    System.out.println(ex.getMessage());
                    setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

                    String msg = ex.getMessage();

                    if (msg.equals("No entry found for ID: null and Class: class uk.ac.ebi.jmzidml.model.mzidml.DBSequence. Make sure the element you are looking for has an ID attribute and is id-mapped!")) {
                        msg = "No dbSequence_ref provided from ProteinDetectionHypothesis, please report this error back to the mzid exporter";
                    }
                    JOptionPane.showMessageDialog(null, msg, "Exception", JOptionPane.ERROR_MESSAGE);
                }

            }
        }.start();

    }

    private void openMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openMenuItemActionPerformed

        progressBarDialog = new ProgressBarDialog(this, true);
        final Thread thread = new Thread(new Runnable() {

            @Override
            public void run() {
                progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                progressBarDialog.setVisible(true);
            }
        }, "ProgressBarDialog");

        int returnVal = fileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            final File mzid_file = fileChooser.getSelectedFile();
            setTitle("ProteoIDViewer   -  " + mzid_file.getPath());
            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    try {
                        if (mzid_file.getPath().endsWith(".gz")) {
                            File outFile = null;
                            FileOutputStream fos = null;

                            GZIPInputStream gin = new GZIPInputStream(new FileInputStream(mzid_file));
                            outFile = new File(mzid_file.getParent(), mzid_file.getName().replaceAll("\\.gz$", ""));
                            fos = new FileOutputStream(outFile);
                            byte[] buf = new byte[100000];
                            int len;
                            while ((len = gin.read(buf)) > 0) {
                                fos.write(buf, 0, len);
                            }
                            fos.close();

                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                            fileName = outFile.getAbsolutePath();
                        } else if (mzid_file.getPath().endsWith(".omx")) {
                            File outFile = null;
                            outFile = new File(fileChooser.getCurrentDirectory(), mzid_file.getName().replaceAll(".omx", ".mzid"));
                            new Omssa2mzid(mzid_file.getPath(), outFile.getPath(), false);
                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                            fileName = outFile.getAbsolutePath();
                        } else if (mzid_file.getPath().endsWith(".xml")) {
                            File outFile = null;
                            outFile = new File(fileChooser.getCurrentDirectory(), mzid_file.getName().replaceAll(".omx", ".mzid"));
                            new Tandem2mzid(mzid_file.getPath(), outFile.getPath());
                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                            fileName = outFile.getAbsolutePath();
                        } else {
                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(mzid_file);
                            fileName = mzid_file.getAbsolutePath();
                        }

                        if (!mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.1.") && !mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.2.")) {
                            progressBarDialog.setVisible(false);
                            progressBarDialog.dispose();
                            JOptionPane.showMessageDialog(null, "The file is not compatible with the Viewer: different mzIdentMl version", "mzIdentMl version", JOptionPane.INFORMATION_MESSAGE);
                            return;
                        }
                        jmzreader = null;
                        createTables();
                        clearSummaryStats();
                        mainTabbedPane.setSelectedIndex(0);
                        secondTab = false;
                        thirdTab = false;
                        fourthTab = false;
                        fifthTab = false;
                        sixthTab = false;
                        loadProteinAmbiguityGroupTable();

                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();

                        String message = "Do you want to load spectrum source file?";

                        int answer = JOptionPane.showConfirmDialog(null, message);
                        if (answer == JOptionPane.YES_OPTION) {
                            JFileChooser fc;
                            //Create a file chooser
                            fc = new JFileChooser();
                            fc.setCurrentDirectory(fileChooser.getCurrentDirectory());

                            fc.addChoosableFileFilter(new SourceFileFilter());
                            int returnVal1 = fc.showOpenDialog(null);

                            if (returnVal1 == JFileChooser.APPROVE_OPTION) {
                                try {
                                    File file = fc.getSelectedFile();
                                    if (file.getAbsolutePath().toLowerCase().endsWith("mgf")) {
                                        jmzreader = new MgfFile(file);
                                        sourceFile = "mgf";
                                        JOptionPane.showMessageDialog(null, file.getName() + " is loaded", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                    } else if (file.getAbsolutePath().toLowerCase().endsWith("mzml")) {
                                        jmzreader = new MzMlWrapper(file);
                                        sourceFile = "mzML";
                                        JOptionPane.showMessageDialog(null, file.getName() + " is loaded", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                    } else {
                                        JOptionPane.showMessageDialog(null, file.getName() + " is not supported", "Spectrum file", JOptionPane.INFORMATION_MESSAGE);
                                    }
                                } catch (JMzReaderException ex) {
                                    System.out.println(ex.getMessage());
                                }
                            }
                        }

                        if (proteinAmbiguityGroupTable.getRowCount() == 0) {
                            JOptionPane.showMessageDialog(null, "There is no protein view for this file", "Protein View", JOptionPane.INFORMATION_MESSAGE);
                        }
                        //loadSummaryStats();

                    } catch (OutOfMemoryError error) {
                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        Runtime.getRuntime().gc();
                        JOptionPane.showMessageDialog(null, "Out of Memory Error.", "Error", JOptionPane.ERROR_MESSAGE);

                        System.exit(0);
                    } catch (Exception ex) {
                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        System.out.println(ex.getMessage());
                        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

                        String msg = ex.getMessage();

                        if (msg.equals("No entry found for ID: null and Class: class uk.ac.ebi.jmzidml.model.mzidml.DBSequence. Make sure the element you are looking for has an ID attribute and is id-mapped!")) {
                            msg = "No dbSequence_ref provided from ProteinDetectionHypothesis, please report this error back to the mzid exporter";
                        }
                        JOptionPane.showMessageDialog(null, msg, "Exception", JOptionPane.ERROR_MESSAGE);
                    }

                }
            }.start();
        }
    }//GEN-LAST:event_openMenuItemActionPerformed
    private void loadDBSequenceTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        while (iterDBSequence.hasNext()) {
            DBSequence dBSequence = iterDBSequence.next();

            String cv = "";
            if (dBSequence.getCvParam() != null) {
                if (dBSequence.getCvParam().size() > 0) {
                    cv = dBSequence.getCvParam().get(0).getValue();
                }
            }

            ((DefaultTableModel) dBSequenceTable.getModel()).addRow(new String[]{
                dBSequence.getId(),
                dBSequence.getAccession(),
                dBSequence.getSeq(),
                cv
            });

        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void loadProtocolData() {
        protocalTextPane.setText("");
        String text = "";
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        HashMap<String, AnalysisSoftware> analysisSoftwareHashMap;
        AnalysisProtocolCollection analysisProtocolCollection;
        List<SpectrumIdentificationProtocol> spectrumIdentificationProtocol;
        ProteinDetectionProtocol proteinDetectionProtocol;

        analysisSoftwareHashMap = new HashMap();

        Iterator<AnalysisSoftware> iterAnalysisSoftware = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.AnalysisSoftware);
        while (iterAnalysisSoftware.hasNext()) {
            AnalysisSoftware analysisSoftware = iterAnalysisSoftware.next();
            analysisSoftwareHashMap.put(analysisSoftware.getId(), analysisSoftware);

        }

        analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);

        spectrumIdentificationProtocol = analysisProtocolCollection.getSpectrumIdentificationProtocol();
        proteinDetectionProtocol = analysisProtocolCollection.getProteinDetectionProtocol();

        if (spectrumIdentificationProtocol != null) {
            text = text + "<html><font color=red>Spectrum Identification Protocol</font><BR>";
            for (int i = 0; i < spectrumIdentificationProtocol.size(); i++) {
                SpectrumIdentificationProtocol spectrumIdentificationProtocol1 = spectrumIdentificationProtocol.get(i);
                Param param = spectrumIdentificationProtocol.get(i).getSearchType();
                if (param != null) {
                    CvParam cVParam = param.getCvParam();
                    if (cVParam != null) {
                        text = text + "<B>Type of search:</B> " + cVParam.getName() + "<BR>";

                    }
                }

                Enzymes enzymes = spectrumIdentificationProtocol.get(i).getEnzymes();
                if (enzymes != null) {
                    List<Enzyme> enzymesList = enzymes.getEnzyme();
                    for (int j = 0; j < enzymesList.size(); j++) {
                        Enzyme enzyme = enzymesList.get(j);
                        if (enzyme != null) {
                            if (enzyme.getEnzymeName() != null) {
                                List<CvParam> cvParamList = enzyme.getEnzymeName().getCvParam();
                                for (int k = 0; k < cvParamList.size(); k++) {
                                    CvParam cvParam = cvParamList.get(k);
                                    text = text + "<B>Enzyme:</B> " + cvParam.getName() + "<BR>";

                                }
                            }
                        }

                    }

                }

                ModificationParams modificationParams = spectrumIdentificationProtocol.get(i).getModificationParams();
                if (modificationParams != null) {
                    List<SearchModification> searchModificationList = modificationParams.getSearchModification();
                    for (int j = 0; j < searchModificationList.size(); j++) {

                        SearchModification searchModification = searchModificationList.get(j);
                        String mod = "";
                        if (searchModification.isFixedMod()) {
                            mod = "(fixed)";
                        } else {
                            mod = "(variable)";
                        }
                        text = text + "<B>Search Modification:</B> Type " + mod + " Residues: " + searchModification.getResidues() + " Mass Delta: " + searchModification.getMassDelta() + "<BR>";

                        if (searchModification.getCvParam() != null) {
                            if (searchModification.getCvParam().get(0) != null) {
                            }
                        }

                    }
                }

                Tolerance tolerance = spectrumIdentificationProtocol.get(i).getFragmentTolerance();
                if (tolerance != null) {
                    String t1 = "";
                    String t2 = "";
                    List<CvParam> cvParamList = tolerance.getCvParam();
                    for (int j = 0; j < cvParamList.size(); j++) {
                        CvParam cvParam = cvParamList.get(j);
                        if (cvParam.getName().equals("search tolerance plus value")) {
                            t1 = t1 + "+" + cvParam.getValue() + " Da ";
                        }
                        if (cvParam.getName().equals("search tolerance minus value")) {
                            t2 = t2 + "-" + cvParam.getValue() + " Da ";
                        }
                    }
                    text = text + "<B>Fragment Tolerance: </B>" + t1 + " / " + t2 + "<BR>";
                }

                Tolerance toleranceParent = spectrumIdentificationProtocol.get(i).getParentTolerance();
                if (toleranceParent != null) {
                    String t1 = "";
                    String t2 = "";
                    List<CvParam> cvParamList = toleranceParent.getCvParam();
                    for (int j = 0; j < cvParamList.size(); j++) {
                        CvParam cvParam = cvParamList.get(j);
                        if (cvParam.getName().equals("search tolerance plus value")) {
                            t1 = t1 + "+" + cvParam.getValue() + " Da ";
                        }
                        if (cvParam.getName().equals("search tolerance minus value")) {
                            t2 = t2 + "-" + cvParam.getValue() + " Da ";
                        }
                    }
                    text = text + "<B>Parent Fragment Tolerance: </B>" + t1 + " / " + t2 + "<BR>";
                }

                String threshold = spectrumIdentificationProtocol.get(i).getThreshold().getCvParam().get(0).getValue();
                if (threshold != null) {
                    text = text + "<B>Threshold: </B>" + threshold + "<BR>";
                }
//                String sw_ref_sip = spectrumIdentificationProtocol1.getAnalysisSoftwareRef();
//                if (sw_ref_sip != null && !sw_ref_sip.equals("")) {
//                      text = text+  "<B>"+"Analysis Software Ref"+": </B> "+      sw_ref_sip+"<BR>";
//                       text = text+  "<B>"+"Analysis Software name"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_sip).getName()+"<BR>";
//                       text = text+  "<B>"+"Analysis Software version"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_sip).getVersion()+"<BR>";
//                       text = text+  "<B>"+"Analysis Software url"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_sip).getUri()+"<BR>";
//                    
//
//                }
            }

        }
        if (proteinDetectionProtocol != null) {
            text = text + "<html><font color=red>Protein Detection Protocol</font><BR>";
            String threshold = proteinDetectionProtocol.getThreshold().getCvParam().get(0).getValue();
            if (threshold != null) {
                text = text + "<B>Threshold: </B>" + threshold + "<BR>";
            }
//           text= text+"<html>Analysis Software Ref<BR>";
//            String sw_ref_pdp = proteinDetectionProtocol.getAnalysisSoftwareRef();
//                        text = text+  "<B>"+"Analysis Software name"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_pdp).getName()+"<BR>";
//                       text = text+  "<B>"+"Analysis Software version"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_pdp).getVersion()+"<BR>";
//                       text = text+  "<B>"+"Analysis Software url"+": </B> "+      analysisSoftwareHashMap.get(sw_ref_pdp).getUri()+"<BR>";
//           

        }
        protocalTextPane.setText(text);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void clearSummaryStats() {

        totalSIRLabelValue.setText("0");
        totalSIIaboveThresholdLabelValue.setText("0");
        totalSIIbelowThresholdLabelValue.setText("0");
        totalSIIaboveThresholdRankOneLabelValue.setText("0");
        percentIdentifiedSpectraLabelValue.setText("0");
        totalPeptidesaboveThresholdLabelValue.setText("0");
        totalPAGsLabelValue.setText("0");
        totalPDHsaboveThresholdLabelValue.setText("0");
        isDecoySiiValue.setText("0");
        isDecoySiiFalseValue.setText("0");

        fpSiiValue.setText("0");
        tpSiiValue.setText("0");
        fdrSiiValue.setText("0");
        siiComboBox.removeAllItems();

        while (fdrPanel.getComponents().length > 0) {
            fdrPanel.remove(0);
        }
        fdrPanel.validate();
        fdrPanel.repaint();

        while (tpEvaluePanel.getComponents().length > 0) {
            tpEvaluePanel.remove(0);
        }
        tpEvaluePanel.validate();
        tpEvaluePanel.repaint();

        while (tpQvaluePanel.getComponents().length > 0) {
            tpQvaluePanel.remove(0);
        }
        tpQvaluePanel.validate();
        tpQvaluePanel.repaint();

    }

    private void loadSummaryStats() {

        Iterator<SpectrumIdentificationList> iterspectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
        while (iterspectrumIdentificationList.hasNext()) {
            SpectrumIdentificationList spectrumIdentificationList1 = iterspectrumIdentificationList.next();

            siiComboBox.addItem(spectrumIdentificationList1.getId());
        }
        loadSpectrumIdentificationList();

    }

    private void loadSpectrumIdentificationList() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        try {
            sIIListPassThreshold.clear();
            sIIListBelowThreshold.clear();
            sIIListPassThresholdRankOne.clear();
            peptideListNonReduntant.clear();
            sIIListIsDecoyFalse.clear();
            sIIListIsDecoyTrue.clear();

            List<SpectrumIdentificationResult> sirListTemp = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationList.class, siiComboBox.getSelectedItem().toString()).getSpectrumIdentificationResult();
            List<SpectrumIdentificationItem> siiListTemp = new ArrayList();
            for (int i = 0; i < sirListTemp.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = sirListTemp.get(i);

                for (int j = 0; j < spectrumIdentificationResult.getSpectrumIdentificationItem().size(); j++) {
                    siiListTemp.add(spectrumIdentificationResult.getSpectrumIdentificationItem().get(j));

                }
            }

            if (siiListTemp.size() > 0) {

                SpectrumIdentificationItem spectrumIdentificationItem = siiListTemp.get(0);
                List<CvParam> cvParamList = spectrumIdentificationItem.getCvParam();
                for (int i = 0; i < cvParamList.size(); i++) {
                    CvParam cvParam = cvParamList.get(i);
                    if (cvParam.getName() != null && !cvParam.getName().equals("")) {
                        boolean isthere = false;
                        for (int j = 0; j < jComboBox1.getItemCount(); j++) {
                            if (jComboBox1.getItemAt(j).equals(cvParam.getName())) {
                                isthere = true;
                                break;
                            }
                        }
                        if (!isthere) {
                            jComboBox1.addItem(cvParam.getName());
                            cvTermHashMap.put(cvParam.getName(), cvParam.getAccession());
                        }

                    }

                }

            }

            for (int j = 0; j < siiListTemp.size(); j++) {
                SpectrumIdentificationItem spectrumIdentificationItem = siiListTemp.get(j);

                if (spectrumIdentificationItem.isPassThreshold()) {
                    sIIListPassThreshold.add(spectrumIdentificationItem);
                    String p_ref = spectrumIdentificationItem.getPeptideRef();

                    // Update FG 18-03-2013 escape XML 
                    Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, StringEscapeUtils.escapeXml(p_ref));

                    if (!peptideListNonReduntant.contains(peptide)) {
                        peptideListNonReduntant.add(peptide);
                    }
                } else {
                    sIIListBelowThreshold.add(spectrumIdentificationItem);
                }

                if (spectrumIdentificationItem.isPassThreshold() && spectrumIdentificationItem.getRank() == 1) {
                    sIIListPassThresholdRankOne.add(spectrumIdentificationItem);
                }

                boolean isdecoy = true;
                List<PeptideEvidenceRef> peRefLst = spectrumIdentificationItem.getPeptideEvidenceRef();
                for (int k = 0; k < peRefLst.size(); k++) {
                    PeptideEvidenceRef peptideEvidenceRef = peRefLst.get(k);
                    PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                    if (!peptideEvidence.isIsDecoy()) {
                        sIIListIsDecoyFalse.add(spectrumIdentificationItem);
                        isdecoy = false;
                        break;
                    }
                }
                if (isdecoy) {
                    sIIListIsDecoyTrue.add(spectrumIdentificationItem);
                }
            }

            int sirSize = mzIdentMLUnmarshaller.getObjectCountForXpath(MzIdentMLElement.SpectrumIdentificationResult.getXpath());
            totalSIRLabelValue.setText("" + sirSize);

            if (sIIListPassThreshold != null) {
                totalSIIaboveThresholdLabelValue.setText(String.valueOf(sIIListPassThreshold.size()));
            }
            if (sIIListBelowThreshold != null) {
                totalSIIbelowThresholdLabelValue.setText(String.valueOf(sIIListBelowThreshold.size()));
            }
            if (sIIListPassThresholdRankOne != null) {
                totalSIIaboveThresholdRankOneLabelValue.setText(String.valueOf(sIIListPassThresholdRankOne.size()));
            }
            if (sIIListPassThresholdRankOne != null && sirSize > 0) {
                double percent = roundTwoDecimals((float) sIIListPassThresholdRankOne.size() * 100 / sirSize);
                percentIdentifiedSpectraLabelValue.setText(String.valueOf(percent) + "%");
            }
            if (peptideListNonReduntant != null) {
                totalPeptidesaboveThresholdLabelValue.setText(String.valueOf(peptideListNonReduntant.size()));
            }
            int pagSize = mzIdentMLUnmarshaller.getObjectCountForXpath(MzIdentMLElement.ProteinAmbiguityGroup.getXpath());

            totalPAGsLabelValue.setText("" + pagSize);

            if (pDHListPassThreshold != null) {
                totalPDHsaboveThresholdLabelValue.setText(String.valueOf(pDHListPassThreshold.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                isDecoySiiValue.setText(String.valueOf(sIIListIsDecoyTrue.size()));
            }
            if (sIIListIsDecoyFalse != null) {
                isDecoySiiFalseValue.setText(String.valueOf(sIIListIsDecoyFalse.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                falsePositiveSii = roundThreeDecimals(sIIListIsDecoyTrue.size());
                fpSiiValue.setText(String.valueOf(falsePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                truePositiveSii = roundThreeDecimals(sIIListPassThreshold.size() - sIIListIsDecoyTrue.size());
                tpSiiValue.setText(String.valueOf(truePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                if (falsePositiveSii + truePositiveSii == 0) {
                    fdrSiiValue.setText("0.0");
                } else {
                    fdrSii = roundThreeDecimals((falsePositiveSii) / (falsePositiveSii + truePositiveSii));
                    fdrSiiValue.setText(String.valueOf(fdrSii));
                }
            }
            //        thresholdValue.setText(spectrumIdentificationProtocol.get(0).getThreshold().getCvParam().get(0).getValue());
        } catch (JAXBException ex) {
            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void fdrPlotActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fdrPlotActionPerformed
    }//GEN-LAST:event_fdrPlotActionPerformed

    private void lookAndFeelDefaultMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeelDefaultMenuItemActionPerformed
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception ex) {
//           
        }

        lookAndFeelDefaultMenuItem.setSelected(true);
        lookAndFeelMetalMenuItem.setSelected(false);
        lookAndFeeMotifMenuItem.setSelected(false);
        repaint();

    }//GEN-LAST:event_lookAndFeelDefaultMenuItemActionPerformed

    private void lookAndFeelMetalMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeelMetalMenuItemActionPerformed
        try {
            UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
        } catch (Exception ex) {
        }
        lookAndFeelDefaultMenuItem.setSelected(false);
        lookAndFeelMetalMenuItem.setSelected(true);
        lookAndFeeMotifMenuItem.setSelected(false);
        repaint();
    }//GEN-LAST:event_lookAndFeelMetalMenuItemActionPerformed

    private void lookAndFeeMotifMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeeMotifMenuItemActionPerformed

        try {
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
        } catch (Exception ex) {
        }
        lookAndFeelDefaultMenuItem.setSelected(false);
        lookAndFeelMetalMenuItem.setSelected(false);
        lookAndFeeMotifMenuItem.setSelected(true);
        repaint();
    }//GEN-LAST:event_lookAndFeeMotifMenuItemActionPerformed

    private void font12ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font12ActionPerformed
        // TODO add your handling code here:
        changeFontSize(12);
        font12.setSelected(true);
        font14.setSelected(false);
        font16.setSelected(false);
    }//GEN-LAST:event_font12ActionPerformed

    private void font14ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font14ActionPerformed
        // TODO add your handling code here:
        changeFontSize(14);
        font12.setSelected(false);
        font14.setSelected(true);
        font16.setSelected(false);
    }//GEN-LAST:event_font14ActionPerformed

    private void font16ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font16ActionPerformed
        // TODO add your handling code here:
        changeFontSize(16);
        font12.setSelected(false);
        font14.setSelected(false);
        font16.setSelected(true);
    }//GEN-LAST:event_font16ActionPerformed

    private void manualCalculateActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_manualCalculateActionPerformed

//            if (!getTitle().endsWith("*")) {
//                setTitle(getTitle() + " *");
//            }
        progressBarDialog = new ProgressBarDialog(this, true);
        final Thread thread = new Thread(new Runnable() {

            @Override
            public void run() {
                progressBarDialog.setTitle("Calculating FDR stats and drawing graphs. Please wait..");
                progressBarDialog.setVisible(true);
            }
        }, "ProgressBarDialog");

        thread.start();

        new Thread("LoadingThread") {

            @Override
            public void run() {
                makeFDRGraphs();
                progressBarDialog.setVisible(false);
                progressBarDialog.dispose();
            }
        }.start();

    }

    public void makeFDRGraphs() {
        try {
            sIIListIsDecoyTrue.clear();
            sIIListIsDecoyFalse.clear();

            List<SpectrumIdentificationResult> sirListTemp = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationList.class, siiComboBox.getSelectedItem().toString()).getSpectrumIdentificationResult();
            List<SpectrumIdentificationItem> siiListTemp = new ArrayList();
            for (int i = 0; i < sirListTemp.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = sirListTemp.get(i);

                for (int j = 0; j < spectrumIdentificationResult.getSpectrumIdentificationItem().size(); j++) {
                    siiListTemp.add(spectrumIdentificationResult.getSpectrumIdentificationItem().get(j));

                }
            }

            boolean isdecoy = true;
            for (int j = 0; j < siiListTemp.size(); j++) {
                SpectrumIdentificationItem spectrumIdentificationItem = siiListTemp.get(j);
                List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();

                for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                    PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);

                    PeptideEvidence peptideEvidence1 = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                    DBSequence dbSeq = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, peptideEvidence1.getDBSequenceRef());

                    if (manualDecoy.isSelected()) {   //Added by ARJ to use value in file if manual decoy is not selected
                        if (!dbSeq.getAccession().startsWith(manualDecoyPrefixValue.getText())) {
                            sIIListIsDecoyFalse.add(spectrumIdentificationItem);
                            isdecoy = false;
                            break;
                        }
                    } else {
                        if (!peptideEvidence1.isIsDecoy()) {
                            sIIListIsDecoyFalse.add(spectrumIdentificationItem);
                            isdecoy = false;
                            break;
                        }
                    }

                }
                if (isdecoy) {
                    sIIListIsDecoyTrue.add(spectrumIdentificationItem);
                }

            }

            if (sIIListIsDecoyTrue != null) {
                isDecoySiiValue.setText(String.valueOf(sIIListIsDecoyTrue.size()));
            }
            if (sIIListIsDecoyFalse != null) {
                isDecoySiiFalseValue.setText(String.valueOf(sIIListIsDecoyFalse.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                falsePositiveSii = roundThreeDecimals(sIIListIsDecoyTrue.size() / Double.valueOf(manualDecoyRatioValue.getText().trim()));
                fpSiiValue.setText(String.valueOf(falsePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                truePositiveSii = roundThreeDecimals(sIIListPassThreshold.size() - sIIListIsDecoyTrue.size());
                tpSiiValue.setText(String.valueOf(truePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                fdrSii = roundThreeDecimals(falsePositiveSii / (falsePositiveSii + truePositiveSii));
                fdrSiiValue.setText(String.valueOf(fdrSii));
            }

            while (fdrPanel.getComponents().length > 0) {
                fdrPanel.remove(0);
            }

            if (manualDecoy.isSelected()) {
                String cvTerm = "";
                boolean order = false;
                if (jComboBox2.getSelectedItem().equals("Better scores are lower")) {
                    order = true;
                }
                cvTerm = cvTermHashMap.get(jComboBox1.getSelectedItem());
                falseDiscoveryRate = new FalseDiscoveryRate(fileName, manualDecoyRatioValue.getText(), manualDecoyPrefixValue.getText(), cvTerm, order);

                falseDiscoveryRate.computeFDRusingJonesMethod();

                // FDR Graph
                XYSeriesCollection datasetFDR = new XYSeriesCollection();
                final XYSeries dataFDR = new XYSeries("FDR", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    if (falseDiscoveryRate.getSorted_evalues().get(i) != 0) {
                        dataFDR.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_estimatedFDR().get(i));
                    }

                }

                final XYSeries dataFDRQvalue = new XYSeries("Q-value", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    if (falseDiscoveryRate.getSorted_evalues().get(i) != 0) {
                        dataFDRQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_qValues().get(i));
                    }

                }

                final XYSeries dataFDRSimple = new XYSeries("Simple FDR", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    if (falseDiscoveryRate.getSorted_evalues().get(i) != 0) {
                        dataFDRSimple.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_simpleFDR().get(i));
                    }

                }

                datasetFDR.addSeries(dataFDR);
                datasetFDR.addSeries(dataFDRQvalue);
                datasetFDR.addSeries(dataFDRSimple);
                final JFreeChart chartFDR = createFDRChart(datasetFDR);
                final ChartPanel chartPanelFDR = new ChartPanel(chartFDR);

                chartPanelFDR.setPreferredSize(new java.awt.Dimension(257, 255));
                fdrPanel.add(chartPanelFDR);
                JScrollPane jFDRPane = new JScrollPane(chartPanelFDR);
                fdrPanel.setLayout(new java.awt.BorderLayout());
                fdrPanel.add(jFDRPane);

                // TP vs Evalue Graph
                XYSeriesCollection datasetTpQvalue = new XYSeriesCollection();
                final XYSeries dataTpQvalue = new XYSeries("TP", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    if (falseDiscoveryRate.getSorted_evalues().get(i) != 0) {
                        dataTpQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getTP().get(i));
                    }

                }

                final XYSeries dataFpQvalue = new XYSeries("FP", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    if (falseDiscoveryRate.getSorted_evalues().get(i) != 0) {
                        dataFpQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getFP().get(i));
                    }

                }

                datasetTpQvalue.addSeries(dataTpQvalue);
                datasetTpQvalue.addSeries(dataFpQvalue);
                final JFreeChart chartTpQvalue = createTpQvalueChart(datasetTpQvalue);
                final ChartPanel chartPanelTpQvalue = new ChartPanel(chartTpQvalue);

                chartPanelTpQvalue.setPreferredSize(new java.awt.Dimension(257, 255));
                tpEvaluePanel.add(chartPanelTpQvalue);
                JScrollPane jTpQvaluePane = new JScrollPane(chartPanelTpQvalue);
                tpEvaluePanel.setLayout(new java.awt.BorderLayout());
                tpEvaluePanel.add(jTpQvaluePane);

                // TP vs Qvalue
                XYSeriesCollection datasetTpQvalueCollection = new XYSeriesCollection();
                final XYSeries dataTpQValueSeries = new XYSeries("data", false);

                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                    dataTpQValueSeries.add(falseDiscoveryRate.getSorted_qValues().get(i), falseDiscoveryRate.getTP().get(i));
                    //System.out.println(falseDiscoveryRate.getSorted_qValues().get(i) + " " + falseDiscoveryRate.getTP().get(i) );
                }

                datasetTpQvalueCollection.addSeries(dataTpQValueSeries);
                final JFreeChart chartTpQvalueChart = createTpQvalue(datasetTpQvalueCollection);
                final ChartPanel chartPanelTpSimpleFDR = new ChartPanel(chartTpQvalueChart);

                chartPanelTpSimpleFDR.setPreferredSize(new java.awt.Dimension(257, 255));
                tpQvaluePanel.add(chartPanelTpSimpleFDR);
                JScrollPane jTpSimpleFDRPane = new JScrollPane(chartPanelTpSimpleFDR);
                tpQvaluePanel.setLayout(new java.awt.BorderLayout());
                tpQvaluePanel.add(jTpSimpleFDRPane);

                repaint();
            }
        } catch (JAXBException ex) {
            Logger.getLogger(ProteoIDViewer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_manualCalculateActionPerformed

    private void manualDecoyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_manualDecoyActionPerformed
        if (manualDecoy.isSelected()) {
            manualDecoyPrefixValue.setEnabled(true);
            manualDecoyRatioValue.setEnabled(true);
            //  manualCalculate.setEnabled(true);
        } else {
            manualDecoyPrefixValue.setEnabled(false);
            manualDecoyRatioValue.setEnabled(false);
            // manualCalculate.setEnabled(false);
        }
}//GEN-LAST:event_manualDecoyActionPerformed

    private void saveMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveMenuItemActionPerformed
//        JFileChooser chooser = new JFileChooser();
//        chooser.setFileFilter(new MzIdentMLFilter());
//
//        chooser.setMultiSelectionEnabled(false);
//        chooser.setDialogTitle("Save");
//
//        File selectedFile;
//
//        int returnVal = chooser.showSaveDialog(this);
//
//        if (returnVal == JFileChooser.APPROVE_OPTION) {
//
//            selectedFile = chooser.getSelectedFile();
//
//            if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
//                selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
//            }
//
//            while (selectedFile.exists()) {
//                int option = JOptionPane.showConfirmDialog(this,
//                        "The  file " + chooser.getSelectedFile().getName()
//                        + " already exists. Replace file?",
//                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);
//
//                if (option == JOptionPane.NO_OPTION) {
//                    chooser = new JFileChooser();
//                    chooser.setFileFilter(new CsvFileFilter());
//                    chooser.setMultiSelectionEnabled(false);
//                    chooser.setDialogTitle("Save");
//
//                    returnVal = chooser.showSaveDialog(this);
//
//                    if (returnVal == JFileChooser.CANCEL_OPTION) {
//                        return;
//                    } else {
//                        selectedFile = chooser.getSelectedFile();
//
//                        if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
//                            selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
//                        }
//                    }
//                } else { // YES option
//                    break;
//                }
//            }
//
//            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
//
//
//            try {
//
//                selectedFile = chooser.getSelectedFile();
//
//                if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
//                    selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
//                }
//
//                if (selectedFile.exists()) {
//                    selectedFile.delete();
//                }
//
//                selectedFile.createNewFile();
//                try {
//                    FileWriter f = new FileWriter(selectedFile);
//                    if (falseDiscoveryRate != null) {
//                        falseDiscoveryRate.writeToMzIdentMLFile(selectedFile.getPath());
//                    }
//                } catch (Exception e) {
//                }
//
//
//
//            } catch (IOException ex) {
//                JOptionPane.showMessageDialog(this,
//                        "An error occured when exporting the spectra file details.",
//                        "Error Exporting",
//                        JOptionPane.ERROR_MESSAGE);
//            }
//
//            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
//        }
    }//GEN-LAST:event_saveMenuItemActionPerformed

    private void siiComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_siiComboBoxActionPerformed
        if (siiComboBox.getSelectedIndex() != -1) {
            loadSpectrumIdentificationList();
        }
    }//GEN-LAST:event_siiComboBoxActionPerformed

    private void exportFDRActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportFDRActionPerformed
        exportFDR();
    }//GEN-LAST:event_exportFDRActionPerformed

    private void mainTabbedPaneMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_mainTabbedPaneMouseClicked
        // compare which tab is selected & boolean

        if (mainTabbedPane.getSelectedIndex() == 1 && !secondTab && mzIdentMLUnmarshaller != null) {
            progressBarDialog = new ProgressBarDialog(this, true);
            final Thread thread = new Thread(new Runnable() {

                @Override
                public void run() {
                    progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                    progressBarDialog.setVisible(true);
                }
            }, "ProgressBarDialog");

            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    loadSpectrumIdentificationResultTable();
                    secondTab = true;
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                }
            }.start();
        }

        if (mainTabbedPane.getSelectedIndex() == 2 && !thirdTab && mzIdentMLUnmarshaller != null) {
            progressBarDialog = new ProgressBarDialog(this, true);
            final Thread thread = new Thread(new Runnable() {

                @Override
                public void run() {
                    progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                    progressBarDialog.setVisible(true);
                }
            }, "ProgressBarDialog");

            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    loadPeptideTable();
                    thirdTab = true;
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                }
            }.start();

        }

        if (mainTabbedPane.getSelectedIndex() == 3 && !fourthTab && mzIdentMLUnmarshaller != null) {
            progressBarDialog = new ProgressBarDialog(this, true);
            final Thread thread = new Thread(new Runnable() {

                @Override
                public void run() {
                    progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                    progressBarDialog.setVisible(true);
                }
            }, "ProgressBarDialog");

            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    loadSpectrumIdentificationResultTable();
                    loadDBSequenceTable();
                    fourthTab = true;
                    progressBarDialog.dispose();
                }
            }.start();

        }

        if (mainTabbedPane.getSelectedIndex() == 4 && !fifthTab && mzIdentMLUnmarshaller != null) {
            progressBarDialog = new ProgressBarDialog(this, true);
            final Thread thread = new Thread(new Runnable() {

                @Override
                public void run() {
                    progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                    progressBarDialog.setVisible(true);
                }
            }, "ProgressBarDialog");

            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    loadSummaryStats();
                    fifthTab = true;
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                }
            }.start();

        }

        if (mainTabbedPane.getSelectedIndex() == 5 && !sixthTab && mzIdentMLUnmarshaller != null) {
            progressBarDialog = new ProgressBarDialog(this, true);
            final Thread thread = new Thread(new Runnable() {

                @Override
                public void run() {
                    progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                    progressBarDialog.setVisible(true);
                }
            }, "ProgressBarDialog");

            thread.start();

            new Thread("LoadingThread") {

                @Override
                public void run() {
                    loadProtocolData();
                    sixthTab = true;
                    progressBarDialog.setVisible(false);
                    progressBarDialog.dispose();
                }
            }.start();

        }
        repaint();

    }//GEN-LAST:event_mainTabbedPaneMouseClicked

    private void psmRankValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_psmRankValueActionPerformed
        loadPeptideTable();
    }//GEN-LAST:event_psmRankValueActionPerformed

    private void jMenuItem1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem1ActionPerformed
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Export Proteins Only");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Export Proteins Only");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                mzidToCsv.useMzIdentMLToCSV(mzIdentMLUnmarshaller, selectedFile.getPath(), "exportProteinsOnly", false);

            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }


    }//GEN-LAST:event_jMenuItem1ActionPerformed

    private void jMenuItem2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem2ActionPerformed
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Export Protein Groups");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Export Protein Groups");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                mzidToCsv.useMzIdentMLToCSV(mzIdentMLUnmarshaller, selectedFile.getPath(), "exportProteinGroups", false);

            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }


    }//GEN-LAST:event_jMenuItem2ActionPerformed

    private void jMenuItem3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem3ActionPerformed
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Export PSMs");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Export PSMs");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                mzidToCsv.useMzIdentMLToCSV(mzIdentMLUnmarshaller, selectedFile.getPath(), "exportPSMs", false);

            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }

    }//GEN-LAST:event_jMenuItem3ActionPerformed

    private void jProteinAmbiguityGroupPanelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jProteinAmbiguityGroupPanelMouseClicked
        // TODO add your handling code here:
    }//GEN-LAST:event_jProteinAmbiguityGroupPanelMouseClicked
    private void changeFontSize(int i) {
        proteinAmbiguityGroupTable.setFont(new Font("Serif", Font.PLAIN, i));
        proteinDetectionHypothesisTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationItemProteinViewTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationResultTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationItemTable.setFont(new Font("Serif", Font.PLAIN, i));
        fragmentationTable.setFont(new Font("Serif", Font.PLAIN, i));
        protocalTextPane.setFont(new Font("Serif", Font.PLAIN, i));
    }

    private JFreeChart createFDRChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "FDR", // chart title
                "log10(e-value)", // x axis label
                "", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
        );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

    private JFreeChart createTpQvalueChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "TP vs FP vs E-value", // chart title
                "log10(e-value)", // x axis label
                "", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
        );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

    private JFreeChart createTpQvalue(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "TP vs Q-value", // chart title
                "Q-value", // x axis label
                "TP value", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
        );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

    private void export(JXTable table, String boxTitle) {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle(boxTitle);

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle(boxTitle);

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                selectedFile.createNewFile();
                try {
                    FileWriter f = new FileWriter(selectedFile);
                    for (int j = 0; j < table.getColumnCount() - 1; j++) {
                        f.write(table.getColumnName(j) + ",");
                    }

                    f.write(table.getColumnName(table.getColumnCount() - 1) + "\r\n");

                    // add the table contents
                    for (int i = 0; i < table.getRowCount(); i++) {
                        for (int j = 0; j < table.getColumnCount() - 1; j++) {
                            f.write(table.getValueAt(i, j) + ",");
                        }

                        f.write(table.getValueAt(i, table.getColumnCount() - 1) + "\r\n");
                    }
                } catch (Exception e) {
                }

            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }

    private void exportFDR() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Export FDR");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Export FDR");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                selectedFile.createNewFile();
                FileWriter f = new FileWriter(selectedFile);
                String outStrHead = "sorted_spectrumResult.get(i)\tsorted_peptideNames.get(i) \t sorted_decoyOrNot.get(i) \t  sorted_evalues.get(i).toString() \t + sorted_scores.get(i).toString() \t estimated_simpleFDR.get(i) \t estimated_qvalue.get(i) \t estimated_fdrscore.get(i) \n";
                f.write(outStrHead);
                for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {

                    String outStr = falseDiscoveryRate.getSorted_spectrumResult().get(i) + "\t"
                            + falseDiscoveryRate.getSorted_peptideNames().get(i) + "\t" + falseDiscoveryRate.getSorted_decoyOrNot().get(i) + "\t"
                            //                    + sorted_evalues.get(i).toString() + "\t" + sorted_scores.get(i).toString() + "\t"
                            + falseDiscoveryRate.getSorted_simpleFDR().get(i) + "\t" + falseDiscoveryRate.getSorted_qValues().get(i) + "\t"
                            + falseDiscoveryRate.getSorted_estimatedFDR().get(i) + "\n";

                    f.write(outStr);
                }

            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {

        if (args != null && args.length == 1) {
            ProteoIDViewerCLI proteoIDViewerCLI = new ProteoIDViewerCLI(args[0]);

        } else {
            java.awt.EventQueue.invokeLater(new Runnable() {

                @Override
                public void run() {

                    new ProteoIDViewer().setVisible(true);
                }
            });
        }
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem aboutMenuItem;
    private javax.swing.JPanel dBSequencePanel;
    private javax.swing.JMenuItem exitMenuItem;
    private javax.swing.JMenuItem exportFDR;
    private javax.swing.JMenu exportMenu;
    private javax.swing.JPopupMenu.Separator exportSeparator2;
    private javax.swing.JPanel fdrPanel;
    private javax.swing.JLabel fdrProteinsLabel;
    private javax.swing.JLabel fdrProteinsValue;
    private javax.swing.JLabel fdrSiiLabel;
    private javax.swing.JLabel fdrSiiValue;
    private javax.swing.JFileChooser fileChooser;
    private javax.swing.JMenu fileMenu;
    private javax.swing.JRadioButtonMenuItem font12;
    private javax.swing.JRadioButtonMenuItem font14;
    private javax.swing.JMenuItem font16;
    private javax.swing.JMenu fontMenu;
    private javax.swing.JLabel fpProteinsLabel;
    private javax.swing.JLabel fpProteinsValue;
    private javax.swing.JLabel fpSiiLabel;
    private javax.swing.JLabel fpSiiValue;
    private javax.swing.JPanel globalStatisticsPanel;
    private javax.swing.JMenu helpMenu;
    private javax.swing.JLabel isDecoySii;
    private javax.swing.JLabel isDecoySiiFalse;
    private javax.swing.JLabel isDecoySiiFalseValue;
    private javax.swing.JLabel isDecoySiiValue;
    private javax.swing.JComboBox jComboBox1;
    private javax.swing.JComboBox jComboBox2;
    private javax.swing.JPanel jExperimentalFilterPanel;
    private javax.swing.JPanel jExperimentalFilterPanel1;
    private javax.swing.JPanel jFragmentationPanel;
    private javax.swing.JPanel jFragmentationPanel1;
    private javax.swing.JPanel jGraph;
    private javax.swing.JPanel jGraph1;
    private javax.swing.JMenuBar jMenuBar;
    private javax.swing.JMenuItem jMenuItem1;
    private javax.swing.JMenuItem jMenuItem2;
    private javax.swing.JMenuItem jMenuItem3;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPeptideEvidencePanel;
    private javax.swing.JPanel jPeptideEvidencePanel1;
    private javax.swing.JPanel jProteinAmbiguityGroupPanel;
    private javax.swing.JEditorPane jProteinDescriptionEditorPane;
    private javax.swing.JPanel jProteinDescriptionPanel;
    private javax.swing.JScrollPane jProteinDescriptionScrollPane;
    private javax.swing.JPanel jProteinDetectionHypothesisPanel;
    private javax.swing.JPanel jProteinInfoPanel;
    private javax.swing.JPanel jProteinSequencePanel;
    private javax.swing.JScrollPane jProteinSequenceScrollPane;
    private javax.swing.JTextPane jProteinSequenceTextPane;
    private javax.swing.JLabel jScientificNameLabel;
    private javax.swing.JLabel jScientificNameValueLabel;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JPopupMenu.Separator jSeparator3;
    private javax.swing.JSeparator jSeparator4;
    private javax.swing.JPanel jSpectrumIdentificationItemPanel;
    private javax.swing.JPanel jSpectrumIdentificationItemPanel1;
    private javax.swing.JPanel jSpectrumIdentificationItemProteinPanel;
    private javax.swing.JPanel jSpectrumIdentificationResultPanel;
    private javax.swing.JPanel jSpectrumPanel;
    private javax.swing.JPanel jSpectrumPanel1;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JSplitPane jSplitPane2;
    private javax.swing.JRadioButtonMenuItem lookAndFeeMotifMenuItem;
    private javax.swing.JRadioButtonMenuItem lookAndFeelDefaultMenuItem;
    private javax.swing.JMenu lookAndFeelMenu;
    private javax.swing.JRadioButtonMenuItem lookAndFeelMetalMenuItem;
    private javax.swing.JPanel mainPanel;
    private javax.swing.JTabbedPane mainTabbedPane;
    private javax.swing.JButton manualCalculate;
    private javax.swing.JCheckBox manualDecoy;
    private javax.swing.JLabel manualDecoyLabel;
    private javax.swing.JLabel manualDecoyPrefix;
    private javax.swing.JTextField manualDecoyPrefixValue;
    private javax.swing.JLabel manualDecoyRatio;
    private javax.swing.JTextField manualDecoyRatioValue;
    private javax.swing.JMenuItem openMenuItem;
    private javax.swing.JMenu optionsMenu;
    private javax.swing.JPanel peptideViewPanel;
    private javax.swing.JLabel percentIdentifiedSpectraLabel;
    private javax.swing.JLabel percentIdentifiedSpectraLabelValue;
    private javax.swing.JPanel proteinDBViewPanel;
    private javax.swing.JPanel proteinViewPanel;
    private javax.swing.JTextPane protocalTextPane;
    private javax.swing.JPanel protocolPanel;
    private javax.swing.JPanel protocolSummaryPanel;
    private javax.swing.JLabel psmRankLabel;
    private javax.swing.JComboBox psmRankValue;
    private javax.swing.JMenuItem saveMenuItem;
    private javax.swing.JComboBox siiComboBox;
    private javax.swing.JLabel siiListLabel;
    private javax.swing.JPanel spectrumViewPanel;
    private javax.swing.JPanel summaryPanel;
    private javax.swing.JLabel totalPAGsLabel;
    private javax.swing.JLabel totalPAGsLabelValue;
    private javax.swing.JLabel totalPDHsLabel;
    private javax.swing.JLabel totalPDHsaboveThresholdLabel;
    private javax.swing.JLabel totalPDHsaboveThresholdLabelValue;
    private javax.swing.JLabel totalPeptidesaboveThresholdLabel;
    private javax.swing.JLabel totalPeptidesaboveThresholdLabelValue;
    private javax.swing.JLabel totalSIILabel;
    private javax.swing.JLabel totalSIIaboveThresholdLabel;
    private javax.swing.JLabel totalSIIaboveThresholdLabelValue;
    private javax.swing.JLabel totalSIIaboveThresholdRankOneLabel;
    private javax.swing.JLabel totalSIIaboveThresholdRankOneLabelValue;
    private javax.swing.JLabel totalSIIbelowThresholdLabel;
    private javax.swing.JLabel totalSIIbelowThresholdLabelValue;
    private javax.swing.JLabel totalSIRLabel;
    private javax.swing.JLabel totalSIRLabelValue;
    private javax.swing.JPanel tpEvaluePanel;
    private javax.swing.JLabel tpProteinsLabel;
    private javax.swing.JLabel tpProteinsValue;
    private javax.swing.JPanel tpQvaluePanel;
    private javax.swing.JLabel tpSiiLabel;
    private javax.swing.JLabel tpSiiValue;
    // End of variables declaration//GEN-END:variables
}
