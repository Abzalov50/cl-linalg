;;; A library on Matrix algebra.
;;; Common Lisp (CL) offers a type `simple-array' with a predicate `arrayp'
;;; to test if an object is an array.
(defpackage :cl-linalg.matrix
  (:nicknames :matrix)
  (:use :cl :cl-linalg.utils)
  (:export :defmatrix :matrix-p
	   :matrix+ :matrix* :matrix** :matrix/ :matrix-map
	   :matrix-loop :matrix-reduce :matrix-augment
	   :matrix-aref :matrix-setf :matrix-flatten
	   :matrix-inverse :gauss-elimination
	   :nullify-coeff :matrix-swap
	   :matrix-max :matrix-min :matrix-extremum
	   :forward-subst :back-subst))

(in-package :cl-linalg.matrix)

(defparameter *padding* 2)
(defparameter *sep* #\Space)

(defun vector= (v1 v2)
  "Return T if the given vectors are equal, otherwise return NIL."
  (every #'true (map 'vector #'equal v1 v2)))

(defun vector-identity-add-p (v)
  "Return T if the vector `v' is an addition identity. Otherwise return NIL."
  (vector= v (identity-add (length v))))

(defun vector-null (v)
  "Return T if all the components of the given vector are 0, otherwise return NIL."
  (vector-identity-add-p v))
#|
(defun vector-identity-mult-p (v)
  "Return T if the vector `v' is a multiplication identity. Otherwise return NIL."
  (vector= v (identity-mult (length v))))
|#

(defun vector-neg (v)
  "Return the opposite vector, in the numerical sense. That is -v"
  (map 'vector #'- v))

(defun vector* (v1 v2)
  "Return the scalar product of vector `v1' by vector `v2'. 
It is assumed that `v1' and `v2' are given as row vectors, and the scalar product is computed by taking the transpose of `v2'. That is scalar-prod = v1 * transpose(v2)."
  (cond ((or (not (vectorp v1)) (not (vectorp v2)))
	 (error "The arguments must be vectors."))
	((not (= (length v1) (length v2)))
	 (error "Vectors must have the same length."))
	((or (vector-null v1) (vector-null v2)) 0)
	((vector-identity-mult-p v1) (reduce #'+ v2))
	((vector-identity-mult-p v2) (reduce #'+ v1))
	(t (reduce #'+ (map 'vector #'* v1 v2)))))

;;; Matrix Algebra
;;; A matrix is a 2D (two-dimensional) array, a vector of vectors.
;;; A (n x m) matrix has `n' rows and `m' columns, or a vector of `n' vectors of length `m'.
;;; `n' is the first dimension and `m' is the second one.
;;; For example:
;;; #((1 2 3) (4 5 6)) is a (2 x 3) array
;;; #((0 -1) (3 3)) is a (2 x 2) array
;;; A vector is just a 1D array.
(defstruct (matrix
	     (:print-function
	      (lambda (struct stream depth)
		(declare (ignore depth))
		(print-matrix stream struct :padding *padding* :sep *sep*)
		(format stream "~&~%[~A x ~A matrix (~A elements)]~%"
			(matrix-nrows struct) (matrix-ncols struct)
			(* (matrix-nrows struct) (matrix-ncols struct))))))
  (dimensions 0)
  (contents #())
  (nrows 0)
  (ncols 0))

;;; Printing functions
(defun type->string (x)
  "Coerce the given data (of any `type') into string.
E.g.: (type->string 125) ==> \"125\"
      (type->string '(1 2 5)) ==> \"(1 2 5)\""
  (typecase x
    (string x)
    (symbol (symbol-name x))
    (t (write-to-string x))))

(defun cols-width (m)
  "Returns a list that contains the width of each column of the given matrix `m'.
The width is the length of the longest element of the column.
To get the length of an element, it is first converted into string."
  (let ((w-lst nil)
	(str-m (matrix-map #'length
			   (matrix-map #'type->string m))))
    (dotimes (i (matrix-ncols m))
      (let ((col (matrix-aref str-m :col i)))
	(setf w-lst
	      (append w-lst
		      (list
		       (if (numberp col)
			   col
			   (matrix-max col)))))))
    w-lst))

(defun print-row (destination row widths &key (padding 1) (sep #\Space))
  "Print the row `row' where each column's width is specified in the list `widths'.
The argument `padding' specify the number of blank space between an element 
and the separator `sep'. Separator can be any character."
  (let ((ncols (if (numberp row)
		   1
		   (matrix-ncols row))))
    (dotimes (i ncols)
      (format destination
	      (concatenate 'string
			   "~"
			   (write-to-string (elt widths i))
			   "<~S~>")
	      (if (= ncols 1)
		  row
		  (aref (matrix-contents row) 0 i)))
      (when (not (= i ncols))
	(format destination
		(concatenate 'string
			     "~" (write-to-string padding)
			     "<~C~>")
		sep)))
    (format destination "~%")))

(defun print-matrix (destination m &key (padding 1) (sep #\Space))
  "Print matrix `m' in tabular format."
  (let ((widths (cols-width m)))
    (dotimes (i (matrix-nrows m))
      (print-row destination (matrix-aref m :row i)
		 widths :padding padding :sep sep))))

(defun copy-array (m)
  "Return a new array `new-m' that contains the same elements as `m'.
That is (eq m new-m) must evaluate to NIL.
By default, (setq new-m m) does not create a new pointer, 
that is (eq m new-m) ==> T."
  (let* ((ndim (array-dimensions m))
	 (new-m (make-array ndim)))
    (cond ((= (length ndim) 1)
	   (dotimes (i (elt ndim 0))
	     (setf (aref new-m i) (aref m i))))
	  ((= (length ndim) 2)
	   (dotimes (i (elt ndim 0))
	     (dotimes (j (elt ndim 1))
	       (setf (aref new-m i j) (aref m i j)))))
	  (t (error "Unsupported dimension of matrix.")))
    new-m))

(defun defmatrix (m &key (inplace nil))
  "Return a matrix structure from the given array (or list) `m'.
When a list or vecteor of length N is given, a matrix of dimension (N 1) is returned,
that is, ALWAYS a column vector."
  (let ((nrows 0) (ncols 0) (contents nil))
    (etypecase m
      (list
       (setq nrows (length m)
	     contents (if inplace m (copy-seq m)))
       (etypecase (car m)
	 (atom
	  (setq ncols 1))
	 (cons
	  (setq ncols (length (car m))))))
      (vector
       (setq nrows (length m)
	     ncols 1
	     contents (if inplace m (copy-seq m))))
      (array
       (setq nrows (array-dimension m 0)
	     ncols (array-dimension m 1)
	     contents (if inplace m (copy-array m))))
      (matrix
       (setq nrows (matrix-nrows m)
	     ncols (matrix-ncols m)
	     contents (if inplace
			  (matrix-contents m)
			  (copy-array (matrix-contents m))))))
    (make-matrix :dimensions (list nrows ncols)
		 :contents contents
		 :nrows nrows
		 :ncols ncols)))

(defun identity-add (size)
  "Return a vector/matrix which is the identity of the vector addition operation.
It is a vector of length `size' filled with 0's."
  (make-array size))

(defun column-vector-p (m)
  "Return T if the given matrix is a column vector, that is a nx1 matrix, otherwise return NIL."
  (and (matrix-p m) (= (matrix-ncols m) 1) (> (matrix-nrows m) 1)))

(defun row-vector-p (m)
  "Return T if the given matrix is a column vector, that is a nx1 matrix, otherwise return NIL."
  (and (matrix-p m) (= (matrix-nrows m) 1) (> (matrix-ncols m) 1)))

(defun dot-product (m1 m2)
  "Return the dot product of the given 1D matrices. `m1' must be a row vector, that is a 1xn matrix, and `m2' must be a column vector, that is a nx1 matrix."
  (cond ((not (row-vector-p m1))
	 (error "The first matrix must be of dimension 1xn."))
	((not (column-vector-p m2))
	 (error "The second matrix must be of dimension nx1."))
	(t (let ((n (matrix-ncols m1)))
	     (if (not (= n (matrix-nrows m2)))
		 (error "The number of columns of the first matrix must be equal to the number of rows of the second one.")
		 (let ((res 0))
		   (dotimes (i n)
		     (incf res
			   (* (matrix-aref m1 :row 0 :col i)
			      (matrix-aref m2 :row i :col 0))))
		   res))))))

(defmacro matrix-loop (m (&key row col) &body body)
  "Loop the given matrix over either the rows, columns or elements (that is over the rows and columns)."
  `(let ((nrows (matrix-nrows ,m))
	 (ncols (matrix-ncols ,m)))
     (cond ((and row col)
	    (dotimes (,row nrows)
	      (dotimes (,col ncols)
		,@body)))
	   (row
	    (dotimes (,row nrows)
	      ,@body))
	   (col
	    (dotimes (,col ncols)
	      ,@body))
	   (t (error "Argument `row' or `col' must be supplied.")))))

(defun matrix-reduce (fn m &key (axis nil))
  "Reduce the given matrix by the given function `fn' (of two arguments) over the given axis. The argument `axis' may be 0, 1 or NIL. If it is 0 (resp. 1) the function `fn' is applied to each two consecutive columns (resp. rows) and the result is a column (resp. row) vector. In this case, `fn' must support 1D matrix arguments. Otherwise, if `axis' is NIL (the default) or any other value, `fn' is applied to each two consecutive elements, and a number is returned. In this case `fn' must support numeric arguments."
  (let ((res nil)
	(nrows (matrix-nrows m))
	(ncols (matrix-ncols m)))
    (cond ((and (numberp axis) (= 0 axis))
	   (setq res (defmatrix (make-array (list nrows 1))))
	   (dotimes (j ncols)
	     (if (zerop j)
		 (setq res (matrix-aref m :col 0))
		 (setq res (funcall fn res (matrix-aref m :col j))))))
	  ((and (numberp axis) (= 1 axis))
	   (setq res (defmatrix (make-array (list 1 ncols))))
	   (dotimes (i nrows)
	     (if (zerop i)
		 (setq res (matrix-aref m :row 0))
		 (setq res (funcall fn res (matrix-aref m :row i))))))
	  (t
	   (dotimes (i nrows)
	     (dotimes (j ncols)
	       (if (and (zerop i) (zerop j))
		   (setq res (matrix-aref m :row 0 :col 0))
		   (setq res
			 (funcall fn res
				  (matrix-aref m :row i :col j))))))))
    res))		 

(defun matrix-map (fn m)
  "Map function `fn' over the elements of matrix `m', 
and return the resulting matrix."
  (let* ((res (make-array (matrix-dimensions m)))
	 (m1 (matrix-contents m))
	 (ncols (matrix-ncols m)))
    (dotimes (i (matrix-nrows m))
      (dotimes (j ncols)
	(setf (aref res i j) (funcall fn (aref m1 i j)))))
    (defmatrix res)))

(defun matrix-map-axis (fn m &optional (row nil) (col nil))
  "Map function `fn' over either the `row'th row or the `col'th column
of matrix `m'. Arguments `row' and `col' exclusive, meaning that only one
should be set. If both are set, an exception is raised."
  (cond ((and row col) (error "Only one of arguments `row' or `col' must be set."))
	(t (let ((res m))
	     (if row
		 (dotimes (i (matrix-ncols m))
		   (setf (aref (matrix-contents res) row i)
			 (funcall fn (aref (matrix-contents m) row i))))
		 (dotimes (i (matrix-nrows m))
		   (setf (aref (matrix-contents res) i col)
			 (funcall fn (aref (matrix-contents m) i col)))))
	     res))))

;;; Operations on matrices
(defun matrix-same-dim-p (m1 m2)
  "Return T if the given matrices have the same dimension,
otherwise return NIL."
  (equal (matrix-dimensions m1) (matrix-dimensions m2)))

(defun matrix-flatten (m)
  "Return a list consisting of the concatenation of the rows of
the given matrix in the right order."
  (let ((nrows (matrix-nrows m))
	 (res nil))
     (dotimes (i nrows)
       (setf res (append res (coerce (matrix-aref m :row i) 'list))))
     res))

(defun matrix-compare (pred &rest args)
  "Apply PREDICATE to every two consecutive matrices of the sequences. 
Return NIL as soon as any invocation of PREDICATE returns NIL, or T if every invocation
is non-NIL.
Two local functions are defined. The first one applies PREDICATE to `only' two matrices.
The second one apply the former recursively to two consecutives matrices."
  (labels ((matrix-comp (pred m1 m2)
	     (if (not (equal (matrix-dimensions m1) (matrix-dimensions m2)))
		 (error "Matrices must have the same dimensions.")
		 (block mapvec
		   (let ((arr1 (matrix-contents m1)) (arr2 (matrix-contents m2)))
		     (dotimes (i (matrix-nrows m1))
		       (dotimes (j (matrix-ncols m1))
			 (when (not (funcall pred (aref arr1 i j) (aref arr2 i j)))
			   (return-from mapvec nil)))))
		   t)))
	   (comp-rec (new-args)
	     (cond ((or (= (length new-args) 1) (null new-args)) t)
		   ((not (funcall #'matrix-comp
				  pred (car new-args) (cadr new-args)))
		    nil)
		   (t (comp-rec (cdr new-args))))))
    (comp-rec args)))

(defun matrix-identity-add-p (m)
  "Return T if the given matrix is an addition identity. Otherwise return NIL."
  (matrix-compare #'= m (defmatrix (identity-add (matrix-dimensions m)))))

(defun matrix-null (m)
  "Return T if all the components of the given matrix are 0, otherwise return NIL."
  (matrix-identity-add-p m))

(defun matrix-square-p (m)
  "Return T if the given MATRIX has the same number of rows and columns,
otherwise return NIL."
  (= (matrix-nrows m) (matrix-ncols m)))

(defun matrix-identity (dimensions)
  "Returns the matrix multiplicative identity of the given dimensions.
It is a matrix where the diagonal elements are 1's and the other elements are 0's.
It only makes sense for square matrices."
  (defmatrix
      (if (not (= (car dimensions) (cadr dimensions)))
	  (error "The number of ROWS and COLUMNS must be equal.")
	  (let ((id-m (make-array dimensions)))
	    (dotimes (i (car dimensions))
	      (setf (aref id-m i i) 1))
	    id-m))))

(defun matrix-identity-p (m)
  "Return T if the given matrix is a multiplication identity. Otherwise return NIL."
  (matrix-compare #'= m (matrix-identity (matrix-dimensions m))))

(defun binary-sub (m1 m2)
  "Return the difference of the two given vectors."
  (cond ((not (= (length m1) (length m2))) (error "Vectors don't have the same size."))
	((matrix-identity-add-p m1) (vector-neg m2))
	((matrix-identity-add-p m2) m1)
	(t (map 'vector #'- m1 m2))))

(defun matrix+- (&rest matrices)
  "Return the sum of the given matrices."
  (labels ((binary-add (m1 m2)
	     ;; Return the sum of the two given matrices
	     (cond ((not (matrix-same-dim-p m1 m2))
		    (error "Matrices must have the same size."))
		   ((matrix-identity-add-p m1) m2)
		   ((matrix-identity-add-p m2) m1)
		   (t (let* ((res-m (defmatrix m1))
			     (arr (matrix-contents res-m))
			     (arr-1 (matrix-contents m1))
			     (arr-2 (matrix-contents m2)))
			(dotimes (i (matrix-nrows m1))
			  (dotimes (j (matrix-ncols m1))
			    (setf (aref arr i j) (+ (aref arr-1 i j) (aref arr-2 i j)))))
			(defmatrix arr))))))
    (if (null matrices)
	0
	(reduce #'binary-add matrices))))

(defun matrix-wise-op (fn default &rest matrices)
  "Apply the given function to the given matrices element-wise."
  (labels ((binary-op (m1 m2)
	     ;; Return the sum of the two given matrices
	     (cond ((not (matrix-same-dim-p m1 m2))
		    (error "Matrices must have the same size."))
		   (t (let* ((res-m (defmatrix m1))
			     (arr (matrix-contents res-m))
			     (arr-1 (matrix-contents m1))
			     (arr-2 (matrix-contents m2)))
			(dotimes (i (matrix-nrows m1))
			  (dotimes (j (matrix-ncols m1))
			    (setf (aref arr i j)
				  (funcall fn
					   (aref arr-1 i j)
					   (aref arr-2 i j)))))
			res-m)))))
    (if (null matrices)
	(if default
	    default
	    (error "Invalid number of arguments"))
	(reduce #'binary-op matrices))))

(defun matrix+ (&rest matrices)
  "Return the sum of the given matrices."
  (apply #'matrix-wise-op #'+ 0 matrices))

(defun matrix/ (&rest matrices)
  "Return the quotient of the given matrices."
  (apply #'matrix-wise-op #'/ nil matrices))

(defun matrix** (&rest matrices)
  "Return the element-wise product of the given matrices."
  (apply #'matrix-wise-op #'* 1 matrices))

(defun matrix-neg (m)
  "Return a matrix where each component is the negative of the corresponding component
in the given matrix, at the same position."
  (cond ((matrix-identity-add-p m) m)
	(t (matrix-map #'- m))))

(defun matrix* (m1 m2)
  "Return the product of the two given matrices."
  (cond ((not (= (matrix-ncols m1) (matrix-nrows m2)))
	 (error "The number of COLUMNS of the first matrix must be equal 
to the number of ROWS of the second one."))
	((or (matrix-identity-add-p m1) (matrix-identity-add-p m2)
	     (matrix-identity-add-p m2)))
	((matrix-identity-p m2) m1)
	((matrix-identity-p m1) m2)
	(t (let* ((nrows (matrix-nrows m1))
		  (ncols (matrix-ncols m2))
		  (res-m (make-array (list nrows ncols))))
	     (dotimes (i nrows)
	       (dotimes (j ncols)
		 (setf (aref res-m i j)
		       (dot-product (matrix-aref m1 :row i)
				    (matrix-aref m2 :col j)))))
	     (defmatrix res-m)))))

(defun matrix-ones (dimensions)
  "Return a matrix in which all the elements are 1's."
  (defmatrix (make-array dimensions :initial-element 1)))

(defun matrix-zeros (dimensions)
  "Return a matrix in which all the elements are 0's."
  (defmatrix (make-array dimensions :initial-element 0)))

(defun matrix-scalar-mult (scalar m)
  "Return a matrix where each element is the product of the given SCALAR
and each element of the given MATRIX."
  (cond ((not (and (matrix-p m) (numberp scalar)))
	      (error "The first argument must be a NUMBER and the second one
must be a MATRIX."))
	(t (matrix-map #'(lambda (x) (* scalar x))
		       m))))

(defun matrix-same-elt-p (eq-pred m)
  "Return T if all the elements of the matrix are equal by the given EQUALITY PREDICATE."
  (let* ((arr (matrix-contents m))
	 (pivot (aref arr 0 0)))
    (if (eq eq-pred #'=)
	(let ((dim (matrix-dimensions m)))
	  (matrix-compare #'= (matrix-ones dim) (matrix-scalar-mult (/ 1 pivot) m)))
	(block comp	
	  (dotimes (i (matrix-nrows m))
	    (dotimes (j (matrix-ncols m))
	      (if (not (funcall eq-pred pivot (aref arr i j)))
		  (return-from comp nil))))
	  t))))

(defun matrix-diag (m)
  "Return the diagonal vector of the given MATRIX.
It makes sense only for square matrices."
  (cond ((not (matrix-square-p m))
	 (error "The number of COLUMNS and ROWS of the given MATRIX must be equal."))
	(t (let ((res-v (make-array (matrix-nrows m)))
		 (arr (matrix-contents m)))
	     (dotimes (i (matrix-nrows m))
	       (setf (aref res-v i) (aref arr i i)))
	     arr))))

(defun matrix-diag-set (diag-v &optional (m nil))
  "Set the given VECTOR as the diagonal of the given MATRIX. If no matrix is given,
return a diagonal matrix.
It makes sense only for square matrices."
  (cond ((and m (not (matrix-square-p m)))
	 (error "The number of COLUMNS and ROWS of the given MATRIX must be equal."))
	((not (null m))
	 (let ((arr (copy-array (matrix-contents m))))
	   (dotimes (i (matrix-nrows m))
	     (setf (aref arr i i) (aref diag-v i)))
	   (defmatrix arr)))
	(t (let* ((len (length diag-v))
		  (arr (identity-add (list len len))))
	     (dotimes (i len)
	       (setf (aref arr i i) (aref diag-v i)))
	     (defmatrix arr)))))

(defun matrix-aref (m &key row col)
  "Return the element, row or column of the given MATRIX at the position 
specified by the ROW and COLUMN arguments."
  (let ((res
	 (cond ((and (not row) (not col))
		m)
	       ((and row col (numberp row) (numberp col))
		(let ((arr (matrix-contents m)))
		  (aref arr row col)))
	       (t (let* ((ncols (matrix-ncols m))
			 (nrows (matrix-nrows m))
			 (row (cond ((not row) (list 0 nrows))
				    ((numberp row) (list row (1+ row)))
				    ((listp row)
				     (if (= 1 (length row))
					 (list (car row) nrows)
					 row))
				    (t (error "Unsupported type"))))
			 (col (cond ((not col) (list 0 ncols))
				    ((numberp col) (list col (1+ col)))
				    ((listp col)
				     (if (= 1 (length col))
					 (list (car col) ncols)
					 col))
				    (t (error "Unsupported type"))))
			 (nc (- (second col) (first col)))
			 (nr (- (second row) (first row)))
			 (new-arr (make-array (list nr nc)))
			 (arr (matrix-contents m)))
		    (do ((i (first row) (1+ i))
			 (p 0 (1+ p)))
			((= i (second row))
			 (defmatrix new-arr))
		      (do ((j (first col) (1+ j))
			   (s 0 (1+ s)))
			  ((= j (second col)))
			(setf (aref new-arr p s) (aref arr i j)))))))))
    (if (and (matrix-p res) (equal (matrix-dimensions res) '(1 1)))
	(let ((arr (matrix-contents res)))
	  (aref arr 0 0))
	res)))

(defun matrix-setf (m newval &key row col)
  "Set the element, row or column of the MATRIX at the given position ROW and COLUMN, to the given VALUE. Even though arguments `row' and `col' are keyword args, it is allowed that only one of them is supplied. The new value `newval' may be a number of a matrix. In the latter case, arguments `row' ad `col' refer to the position or index of the element of the matrix `m' that must be replaced by the element of the matrix `newval' at the position row=0, col=0.
m: matrix, newval: number | matrix, row: integer, col: integer. "
  (cond ((and (not row) (not col))
	 m)
	((not (and row col))
	 (error "If argument ROW is supplied, argument COL must be also be supplied, vice-versa."))
	(t
	 (let ((arr (matrix-contents m)))
	   (cond ((numberp newval)
		  (setf (aref arr row col) newval))
		 ((matrix-p newval)
		  (let* ((ncols (matrix-ncols newval))
			 (nrows (matrix-nrows newval))
			 (newarr (matrix-contents newval)))
		    (dotimes (i nrows)
		      (dotimes (j ncols)
			(setf (aref arr (+ row i) (+ col j))
			      (aref newarr i j))))))
		 (t (error "The new value must be either a NUMBER or a MATRIX.")))
	   m))))

(defun matrix-T (m)
  "Return the transpose of the given matrix.
It makes sense only for square matrices."
  (cond ((and (not (matrix-vector-p m))
	      (not (matrix-square-p m)))
	 (error "The number of COLUMNS and ROWS of the given MATRIX must be equal or the matrix must be a 1D matrix"))
	((and (not (matrix-vector-p m))
	      (or (matrix-same-elt-p #'equalp m)
		  (matrix-identity-p m)))
	 (defmatrix m))
	((and (not (matrix-vector-p m))
	      (matrix-identity-add-p m))
	 (defmatrix m))
	(t
	 (let* ((new-m (defmatrix m))
		(arr (matrix-contents new-m))
		(nrows (matrix-nrows new-m))
		(ncols (matrix-ncols new-m)))
	   (cond ((matrix-column-p m)
		  (update-dim new-m (list 1 nrows))
		  (setf (matrix-contents new-m)
			(let ((arr-1
			       (make-array (matrix-dimensions new-m))))
			  (dotimes (j nrows)
			    (setf (aref arr-1 0 j)
				  (aref arr j 0)))
			  arr-1)))
		 ((matrix-row-p m)
		  (update-dim new-m (list ncols 1))
		  (setf (matrix-contents new-m)
			(let ((arr-1
			       (make-array (matrix-dimensions new-m))))
			  (dotimes (i ncols)
			    (setf (aref arr-1 i 0)
				  (aref arr 0 i)))
			  arr-1)))
		 (t (dotimes (i nrows)
		      (dotimes (j ncols)
			;; leave the diagonal unchanged
			(when (not (= i j))  
			  (setf (aref arr i j) (aref arr j i)))))))
	   new-m))))

(defun matrix-column-p (m)
  "Return T is the givem matrix is a column array, that is, its second dimension is 1,
otherwise return NIL."
  (= (matrix-ncols m) 1))

(defun matrix-row-p (m)
  "Return T is the givem matrix is a row array, that is, its first dimension is 1,
otherwise return NIL."
  (= (matrix-nrows m) 1))

(defun matrix-vector-p (m)
  "Return T is the givem matrix is a vector, that is, one of its dimensions is 1,
otherwise return NIL."
  (or (= (matrix-nrows m) 1) (= (matrix-ncols m) 1)))

(defun matrix-symetric-p (m)
  "Return T if the given MATRIX is symmetric, otherwise return NIL.
It makes sense only for square matrices."
  (cond ((not (matrix-square-p m))
	 (error "The number of COLUMNS and ROWS of the given MATRIX must be equal."))
	(t (matrix-compare #'equalp m (matrix-T m)))))

(defun update-dim (m dim)
  (setf (matrix-nrows m) (car dim)
	(matrix-ncols m) (cadr dim)
	(matrix-dimensions m) dim))

(defun matrix-augment (m1 m2 &key (axis 0))
  "Augment matrix `m1' with matrix `m2' along the given axis defaulted to 0."
  (when (or (and (= axis 1)
		 (not (= (matrix-ncols m1) (matrix-ncols m2))))
	    (and (= axis 0)
		 (not (= (matrix-nrows m1) (matrix-nrows m2)))))
      (error "The two given matrices must have conforming dimensions."))
  (let* ((nrows (matrix-nrows m1)) (ncols (matrix-ncols m1))
	 (nrows-2 (matrix-nrows m2)) (ncols-2 (matrix-ncols m2)))
    (cond ((= axis 0)
	   (let* ((new-m (matrix-zeros (list nrows (+ ncols ncols-2)))))
	     (dotimes (j ncols)
	       (matrix-setf new-m
			    (matrix-aref m1 :col j) :row 0 :col j))
	     (dotimes (j ncols-2)
	       (matrix-setf new-m
			    (matrix-aref m2 :col j)
			    :row 0
			    :col (+ ncols j)))
	     (update-dim new-m (list nrows (+ ncols ncols-2)))
	     new-m))
	  ((= axis 1)
	   (let* ((new-m (matrix-zeros (list (+ nrows nrows-2) ncols))))
	     (dotimes (i nrows)
	       (matrix-setf new-m (matrix-aref m1 :col j)
			    :row i :col 0))
	     (dotimes (i nrows-2)
	       (matrix-setf new-m
			    (matrix-aref m2 :row i)
			    :row (+ nrows i) :col 0))
	     (update-dim new-m (list (+ nrows nrows-2) ncols))
	     new-m))
	  (t (error "Matrix does not have this axis.")))))
 
(defun matrix-extremum (m &key kind (key #'identity) (test #'identity))
  "Return the extremum value (and its index) of the given matrix."
  (let* ((arr (matrix-contents m))
	 (extremum nil)
	 (row nil)
	 (col nil)
	 (fn (if (eq kind 'max)
		 #'>
		 #'<))
	 (flag nil))
    (dotimes (i (matrix-nrows m))
      (dotimes (j (matrix-ncols m))
	(let ((val (aref arr i j)))
	  (when (funcall test val)
	    (when (null extremum)
	      (setf extremum val
		    row i
		    col j))
	    (when (funcall fn (funcall key val) extremum)
	      (setf extremum val
		    row i
		    col j))))))
    (values extremum
	    (cons row col))))

(defun matrix-max (m &key (key #'identity) (test #'identity))
  "Return the maximum value (and its index) of the given matrix."
  (matrix-extremum m :kind 'max :key key :test test))

(defun matrix-min (m &key (key #'identity) (test #'identity))
  "Return the minimum value (and its index) of the given matrix."
  (matrix-extremum m :kind 'min :key key :test test))

(defun matrix-swap (m idx-1 idx-2 &key (axis 0))
  "Swap vectors of the given indices along the given axis.
If `inplace' is set to T, the given is DESTRUCTIVELY modified and
returned, otherwise a new matrix is created."
  (cond ((= axis 0)
	 (let ((temp (matrix-aref m :row idx-1)))
	   (matrix-setf m (matrix-aref m :row idx-2) :row idx-1)
	   (matrix-setf m temp :row idx-2)))
	((= axis 1)
	 (let ((temp (matrix-aref m :col idx-1)))
	   (matrix-setf m (matrix-aref m :col idx-2) :col idx-1)
	   (matrix-setf m temp :col idx-2)))
	(t (error "Matrix does not have this axis."))))

(defun gauss-elimination (A b)
  "Perform a Gauss elimination on the given matrix with the second argument being the right-hand side of the system of linear equations of the form `Ax=b'"
  (let ((nrows (matrix-nrows A))
	(ncols (matrix-ncols A))
	(tol 1e-6)
	(A (matrix-augment A b :axis 0)))
    (labels ((select-pivot (A col)
	       (do* ((i col (1+ i))
		     (flag nil (= i nrows)))
		    ((or (= i nrows)
			 (> (abs (matrix-aref A :row i :col col))
			    tol))
		     (values (if flag
				 nil
				 (matrix-aref A :row i :col col))
			     i flag)))))
      (block outer
	(dotimes (j ncols)
	  (block inner
	    (when (= (1+ j) nrows)
	      (return-from outer A))
	    (multiple-value-bind (pivot p-idx flag)
		(select-pivot A j)
	      (when flag (return-from inner))
	      (setf A (matrix-swap A j p-idx))
	      (do* ((i (1+ j) (1+ i)))
		   ((= i nrows))
		(let ((row-modif (matrix-aref A :row i))
		      (lead (matrix-aref A :row i :col j)))
		  (tagbody
		     (when (zerop lead)
		       (go continue))
		     (setf A (matrix-setf
			      A
			      (matrix+ row-modif
				       (matrix-neg
					(matrix-scalar-mult
					 (/ lead pivot)
					 (matrix-aref A :row p-idx))))
			      :row i))
		   continue))))))))))

;;(setf A (defmatrix (make-array '(3 3)
;;			       :initial-contents '((1 2 3) (-3 -2 -1) ;;(4 4 4)))))

;;(setf b (defmatrix (make-array '(3 1) :initial-contents '((1) (2) (3)))))
#|
(defun gauss-seidel (A b &optional x)
  (let ((A (matrix-augment A b))
	(x (if x
	       x
	       (make-array (matrix-dimensions b)
			   :initial-element 1))))
    (labels ((rec (A x)
	       (cond ((convergedp x)
		      x)
		     (t
		      (dotimes (i (matrix-nrows b))
			(let ((imprv (matrix-aref b :row i)))
			  (dotimes (j (matrix-ncols A))
			    (decf imprv
				  (if (= i j)
				      0
				      (* (matrix-aref A :row i :col j)
					 (aref x i 0))))))
			  (setf (aref x i)
				(/ imprv
				   (aref x i)))))))))))
|#
(defun select-pivot (m col)
  ;; find the row of the maximum absolute coefficient
  ;; of the given `col'
  (let ((p-idx col)
	(p-val (matrix-aref m :row col :col col))
	(nrows (matrix-nrows m)))
    (do* ((i col (1+ i)))
	 ((= i nrows)
	  (values p-val p-idx))
      (let ((curr (matrix-aref m :row i :col col)))
	(when (> (abs curr) (abs p-val))
	  (setf p-idx i
		p-val curr))))))

(defun matrix-inverse (m)
  "Compute matrix inverse using Gauss-Jordan elimination on the given matrix."
  (let* ((nrows (matrix-nrows m))
	 (ncols (matrix-ncols m))
	 (tol 1e-6)
	 (I (matrix-identity (matrix-dimensions m)))
	 (m (matrix-augment m I :axis 0)))
    (block outer
      (dotimes (j ncols)
	(block inner
	  (when (= j nrows)
	    (return-from outer m))
	  (multiple-value-bind (pivot p-idx)
	      (select-pivot m j)
	    (when (< (abs pivot) tol)
	      (error "The given matrix is not invertible!!"))
	    ;; Interchange rows
	    (setf m (matrix-swap m j p-idx))
	    ;; Nullify other coefficients of current column
	    (setf m (nullify-coeff m p-idx j))))))
    (matrix-aref m :col (list ncols (* 2 ncols)))))

(defun nullify-coeff (m row col)
  "Return matrix in which all the coefficients, except the pivot, are nullified using elementary row operations."
  (let* ((pivot-row (matrix-aref m :row row))
	 (nrows (matrix-nrows m))
	 (pivot (matrix-aref m :row row :col col)))
    ;; Normalize pivotal row so that pivot = 1
    (matrix-setf m
		 (matrix-scalar-mult (/ 1 pivot)
				     pivot-row)
		 :row row :col 0)
    ;;(print m)
    (setf pivot-row (matrix-aref m :row row))
    ;;(format t "~&pivot-row:~%~A~%" pivot-row)
    (do* ((i 0 (1+ i)))
	 ((= i nrows) m)
      (when (not (= i row))
	(let ((row-modif (matrix-aref m :row i))
	      (lead (matrix-aref m :row i :col col)))
	  (tagbody
	     (when (zerop lead)
	       (go continue))
	     (matrix-setf
	      m
	      (matrix+ row-modif
		       (matrix-neg (matrix-scalar-mult
				    lead pivot-row)))
	      :row i :col 0)
	     ;;(print m)
	   continue))))))

(defun back-subst (m)
  "Given a (m x m+1) matrix, where the block (m x m) is upper triangular and the remaining vector (m x 1) is the right-hand side of a linear equation system `Ax = b', return the solution `x' using back substitution."
  (let* ((nrows (1- (matrix-nrows m)))
	 (ncols (1- (matrix-ncols m)))
	 (x (defmatrix (make-array (list (1+ nrows) 1))))
	 (x-last (/ (matrix-aref m :row nrows
				   :col ncols)
		      (matrix-aref m :row nrows
				   :col (- ncols 1)))))
    (matrix-setf x x-last
		 :row nrows :col 0)
    (do ((i (- nrows 1) (decf i)))
	((< i 0) x)
      (let ((rhs (matrix-aref m :row i :col ncols))
	    (pivot (matrix-aref m :row i :col i))
	    (coeff-x-known
	     (matrix-aref m :row i :col (list (1+ i) ncols)))
	    (x-known (matrix-aref x :row (list (1+ i) (1+ nrows)))))
	(matrix-setf
	 x
	 (/ (- rhs
	       (if (numberp x-known)
		   (* coeff-x-known x-known)
		   (dot-product coeff-x-known x-known)))
	    pivot)
	 :row i :col 0)))
    x))

(defun forward-subst (m)
  "Given a (m x m+1) matrix, where the block (m x m) is lower triangular and the remaining vector (m x 1) is the right-hand side of a linear equation system `Ax = b', return the solution `x' using forward substitution."
  (let* ((nrows (1- (matrix-nrows m)))
	 (ncols (1- (matrix-ncols m)))
	 (x (defmatrix (make-array (list (1+ nrows) 1))))
	 (x-first (/ (matrix-aref m :row 0 :col ncols)
		      (matrix-aref m :row 0 :col 0))))
    (matrix-setf x x-first
		 :row 0 :col 0)
    (do ((i 1 (incf i)))
	((> i nrows) x)
      (let ((rhs (matrix-aref m :row i :col ncols))
	    (pivot (matrix-aref m :row i :col i))
	    (coeff-x-known
	     (matrix-aref m :row i :col (list 0 i)))
	    (x-known (matrix-aref x :row (list 0 i))))
	(matrix-setf
	 x
	 (/ (- rhs
	       (if (numberp x-known)
		   (* coeff-x-known x-known)
		   (dot-product coeff-x-known x-known)))
	    pivot)
	 :row i :col 0)))
    x))
