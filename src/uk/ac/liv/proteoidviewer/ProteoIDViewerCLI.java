
package uk.ac.liv.proteoidviewer;

/**
 *
 *  @author Fawaz Ghali
 */
public class ProteoIDViewerCLI {
    private  String mzidFile;
    private  ProteoIDViewer proteoIDViewer=null;
   
    public ProteoIDViewerCLI(final String mzidFile){
        this.mzidFile=mzidFile;
        java.awt.EventQueue.invokeLater(new Runnable() {

            @Override
            public void run() {

               proteoIDViewer=  new ProteoIDViewer();
               proteoIDViewer.setVisible(true);
               proteoIDViewer.openCLI(mzidFile);
               
               
            }
        });
    }
    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        ProteoIDViewerCLI proteoIDViewerCLI = new ProteoIDViewerCLI(args[0]);
        
        
    } 
}
