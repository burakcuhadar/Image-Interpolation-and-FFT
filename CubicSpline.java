import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 *
 * @author braeckle
 *
 */
public class CubicSpline implements InterpolationMethod {

    /** linke und rechte Intervallgrenze x[0] bzw. x[n] */
    double a, b;

    /** Anzahl an Intervallen */
    int n;

    /** Intervallbreite */
    double h;

    /** Stuetzwerte an den aequidistanten Stuetzstellen */
    double[] y;

    /** zu berechnende Ableitunge an den Stuetzstellen */
    double yprime[];

    /**
     * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
     * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
     * die Ableitungen an den Stellen x[0] und x[n] = 0.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        this.a = a;
        this.b = b;
        this.n = n;
        h = ((double) b - a) / (n);

        this.y = Arrays.copyOf(y, n + 1);

        /* Randbedingungen setzten */
        yprime = new double[n + 1];
        yprime[0] = 0;
        yprime[n] = 0;

        /* Ableitungen berechnen. Nur noetig, wenn n > 1 */
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * getDerivatives gibt die Ableitungen yprime zurueck
     */
    public double[] getDerivatives() {
        return yprime;
    }

    /**
     * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
     * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
     */
    public void setBoundaryConditions(double yprime0, double yprimen) {
        yprime[0] = yprime0;
        yprime[n] = yprimen;
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
     * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
     * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
     * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
     * Membervariable yprime gespeichert.
     *
     * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
     * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
     * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
     * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
     * berechnet werden muessen.
     */
    public void computeDerivatives() {       
        
        if( n == 2) {
        	yprime[1] = (3.0 / h * (y[2] - y[0] - h / 3.0 * yprime[0]) - yprime[2]) / 4;
        	return;
        }
        
        // Thomas-Algorithmus
        double[] cprime = new double[ n - 2 ];
        double[] dprime = new double[ n - 1];
        
        cprime[0] = 1.0 / 4.0;
        for(int i=1; i<cprime.length; i++)
        	cprime[i] = 1.0 / (4.0 - cprime[i-1]);
        
        
        dprime[0] = 3.0 / h * (y[2] - y[0] - h / 3.0 * yprime[0]) / 4.0;
        for(int i=1; i<dprime.length-1; i++)
        	dprime[i] = (3.0 / h * (y[i+2] - y[i]) - dprime[i-1]) / (4.0 - cprime[i-1]);
              
        dprime[n - 2] = (3.0 / h * (y[n] - y[n-2] - h / 3.0 * yprime[n]) - dprime[n-3]) / (4.0 - cprime[n-3]);
    
        yprime[n - 1] = dprime[n - 2];
        for(int i=n-2; i>=1; i--)
        	yprime[i] = dprime[i-1] - cprime[i-1] * yprime[i+1]; 
    
    }

    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
     * und das entsprechende kubische Hermite-Polynom ausgewertet.
     */
    @Override
    public double evaluate(double z) {
        if(z < a)
        	return y[0];
        
        if(z > b)
        	return y[n];
        
        // Find the interval
        double x_i = a;
        while(x_i < z) {
        	x_i += h;
        }
        if(x_i != a)
        	x_i -= h;
                
        int i = (int) ((x_i - a) / h);
        
        double t = (z - x_i) / h;
        
        return	y[i] 			* (1 - 3*t*t + 2*t*t*t)
        	 +  y[i+1] 			* (3*t*t - 2*t*t*t)
        	 +  h * yprime[i] 	* (t - 2*t*t + t*t*t)
        	 +  h * yprime[i+1]	* (-t*t + t*t*t);
    }	
    
}
