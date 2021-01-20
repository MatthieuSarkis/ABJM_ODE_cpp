#pragma once

template <class T>
class MathFunction
{
public:
  MathFunction();
  virtual ~MathFunction();
  virtual void operator()(double t, T q, T qp) = 0;
};

template <class T>
MathFunction<T>::MathFunction()
{

}

template <class T>
MathFunction<T>::~MathFunction()
{

}
