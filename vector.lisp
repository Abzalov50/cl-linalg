(defun vector-compare (pred &rest args)
  "Apply PREDICATE to every two consecutive vectors of the sequences. 
Return NIL as soon as any invocation of PREDICATE returns NIL, or T if every invocation
is non-NIL.
Two local functions are defined. The first one applies PREDICATE to `only' two vectors.
The second one apply the former recursively to two consecutives vectors."
  (labels ((vec-comp (pred v1 v2)
	     (if (not (equal (length v1) (length v2)))
		 (error "Vectors must have the same size.")
		 (block mapvec
		   (map 'vector #'(lambda (elt1 elt2)
				    (when (not (funcall pred elt1 elt2))
				      (return-from mapvec nil)))
			v1 v2)
		   t)))
	   (comp-rec (new-args)
	     (cond ((or (= (length new-args) 1) (null new-args)) t)
		   ((not (funcall #'vec-comp
				  pred (car new-args) (cadr new-args)))
		    nil)
		   (t (comp-rec (cdr new-args))))))
    (comp-rec args)))

(defun vector- (&rest vectors)
  "Return the substraction of the vectors given as arguments."
  (labels ((binary-sub (v1 v2)
	     ;; Return the difference of the two given vectors
	     (cond ((not (= (length v1) (length v2)))
		    (error "Vectors must have the same size."))
		   ((vector-identity-add-p v1) (vector-neg v2))
		   ((vector-identity-add-p v2) v1)
		   (t (map 'vector #'- v1 v2)))))
    (if (null vectors)
	0
	(reduce #'binary-sub vectors))))

(defun vector-scalar-mult (a v)
  "Multiply the vector `v' by the scalar `a', and return the product."
  (cond ((or (not (numberp a)) (not (vectorp v))) (error "Arguments are not of the correct type."))
	((zerop a) (identity-add (length v)))
	((= a 1) v)
	((vector-null v) v)
	(t (map 'vector #'(lambda (x) (* a x)) v))))

(defun vector-cross-product (v1 v2)
  "Return the cross product of the given vectors, `v1' and `v2'.
It is assumed that the two vectors are given as row vectors, and the cross product is computed
by taking the transpose of `v1'. That is scalar-prod = transpose(v1) * v2.
This gives a tabular data (or matrix or 2D-array) as a result."
  (let ((n (length v1)))
    (cond ((or (not (vectorp v1)) (not (vectorp v2)))
	   (error "The arguments must be vectors."))
	  ((not (= n (length v2))) (error "Vectors must have the same length."))
	  ((or (vector-null v1) (vector-null v2))
	   (make-array (list n n) :initial-element 0))
	  (t (let ((res (make-array (list n n))))
	       (dotimes (i n)
		 (dotimes (j n)
		   (setf (aref res i j) (* (aref v1 i) (aref v2 j)))))
	       res)))))

(defun vector-norm (v)
  "Return the Euclidean norm (or module) of the given vector."
  (if (not (vectorp v))
      (error "The argument must be a vector.")
      (sqrt (reduce #'+ (map 'vector #'square v)))))

(defun vector-distance (v1 v2)
  "Return the distance the distance between two vectors in the Euclidean space."
  (if (or (not (vectorp v1)) (not (vectorp v2)))
      (error "The arguments must be vectors.")
      (vector-norm (vector- v1 v2))))

;;; Statistics functions
(defun vector-avg (v)
  "Return the aithmetical average of the number in the vector.
The vector is viewed as a plain set."
  (if (not (vectorp v))
      (error "The argument must be a vector.")
      (float (/ (reduce #'+ v) (length v)))))

(defun vector-median (v)
  "Return the median of the vector. It is the element that split the vector (set)
into two subsets of the same size."
  (if (not (vectorp v))
      (error "The argument must be a vector.")
      (let* ((n (length v))
	     (mean-idx (/ n 2))
	     (w (sort v #'<))
	     (median (if (evenp n)
			 (/ (+ (aref w (1- mean-idx)) (aref w mean-idx)) 2.0)
			 (aref v (floor n 2)))))
	median)))

(defun vector-min (v)
  "Return the minimum value of the vector (set) given in argument."
  (apply #'min (coerce v 'list)))

(defun vector-max (v)
  "Return the maximum value of the vector (set) given in argument."
  (apply #'max (coerce v 'list)))

(defun vector-stddev (v)
  "Return the standard deviation of the set `v' assuming that it follows a uniform distribution."
  (let ((avg (vector-avg v))
	(n (length v)))
    (sqrt (/ (reduce #'+ (map 'vector
			      #'(lambda (x) (square (- x avg)))
			      v))
	     n))))
