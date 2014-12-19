package cgl.imr.samples.dacidr.wdasmacof.vary;

/**
 * @author Yang Ruan, yangruan@cs.indiana.edu
 * 
 *         Exception for matrix operations.
 * 
 */
public class MatrixException extends Exception {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public MatrixException() {
		super();
	}

	public MatrixException(String message) {
		super(message);
	}

	public MatrixException(String message, Throwable cause) {
		super(message, cause);
	}

	public MatrixException(Throwable cause) {
		super(cause);
	}
}
