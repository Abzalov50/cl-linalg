(defpackage :cl-linalg.utils
  (:use :cl)
  (:export :true :square))

(in-package :cl-linalg.utils)

(defun true (arg)
  (if arg
      t
      nil))

(defun square (x)
  (* x x))
