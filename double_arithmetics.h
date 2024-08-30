/**
\file
\brief Provides means of comparing floating point values concerning deviation due to calculating errors
*/

#ifndef DOUBLE_ARITHMETICS_H
#define DOUBLE_ARITHMETICS_H

/// Maximum difference between floating point values for them to be considered equal
const double DOUBLE_DEVIATION = 1e-20;

// TODO: add and fix docs
/**
 @brief Defines whether a floating point is zero
 \param[in] double \b a

 \return \b true if a=0
 \return \b false otherwise
 */
bool is_zero(const double a);

/// Defines whether floating point values are equal
/**
    \param[in] double <b>a, b</b>

    \return \b true if a=b
    \return \b false otherwise
*/
bool are_equal(const double a, const double b);

double min(const double a, const double b);

double max(const double a, const double b);


#endif // DOUBLE_ARITHMETICS_H
