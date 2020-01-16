;;; -*- mode: lisp -*-
(defpackage :cl-linalg-system
  (:use :cl :asdf))
(in-package :cl-linalg-system)

(defsystem cl-linalg
  :author "Arnold N'GORAN"
  :licence "LLGPL"
  :components ((:file "utils")
	       (:file "matrix"
		      :depends-on ("utils")))
  :depends-on ())
