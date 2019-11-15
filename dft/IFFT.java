package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
    	int n = c.length;
    	
    	Complex[] v = new Complex[n];
    	
    	if(n == 1) {
    		v[0] = c[0];
    		return v;
    	}
    	
    	int m = n/2;
    	
    	Complex[] c_even = new Complex[n/2];
    	Complex[] c_odd = new Complex[n/2];
    	
    	for(int i=0; i<n; i++)
    		if(i % 2 == 0)
    			c_even[i/2] = c[i];
    		else
    			c_odd[(i-1)/2] = c[i];
    	
    	Complex[] z1 = ifft(c_even);
    	Complex[] z2 = ifft(c_odd);    	
    	
    	Complex omega = Complex.fromPolar(1, 2 * Math.PI / n);
    	
    	for(int j=0; j<m; j++) {
    		v[j] = z1[j].add( omega.power(j).mul( z2[j] ));
    		v[m+j] = z1[j].sub( omega.power(j).mul(z2[j]));
    	}
    	
    	return v;
    }
        
}
