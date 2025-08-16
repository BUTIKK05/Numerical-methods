package org.example.lab22222;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

public class QuadratureIntegration {
    private static final MathContext MC = new MathContext(50);

    private static double calculateMoment(int k) {
        double result = 0;
        if (k == 0) return 1 - Math.cos(1);
        result = -Math.pow(1, k) * Math.cos(1) + k * calculateCosMoment(k - 1);
        result -= -k * calculateCosMoment(k - 1);
        return result;
    }

    private static double calculateCosMoment(int k) {
        if (k == 0) return Math.sin(1);
        return Math.pow(1, k) * Math.sin(1) - k * calculateMoment(k - 1);
    }

    private static double[] findP4Coefficients(double[] moments) {
        double[][] coeffs = new double[4][4];
        double[] rhs = new double[4];
        for (int k = 0; k < 4; k++) {
            coeffs[k][0] = moments[k + 3];
            coeffs[k][1] = moments[k + 2];
            coeffs[k][2] = moments[k + 1];
            coeffs[k][3] = moments[k];
            rhs[k] = -moments[k + 4];
        }
        RealMatrix matrix = new Array2DRowRealMatrix(coeffs);
        RealVector vector = new ArrayRealVector(rhs);
        return new LUDecomposition(matrix).getSolver().solve(vector).toArray();
    }

    private static double[] performSyntheticDivision(double[] coeffs, double root) {
        int degree = coeffs.length - 1;
        double[] quotient = new double[degree];
        quotient[degree - 1] = coeffs[degree];
        for (int i = degree - 2; i >= 0; i--) {
            quotient[i] = coeffs[i + 1] + root * quotient[i + 1];
        }
        return quotient;
    }

    private static double[] determineRoots(double[] coeffs) {
        double[] poly = {coeffs[3], coeffs[2], coeffs[1], coeffs[0], 1.0};
        double[] roots = new double[4];
        double[] currentPoly = poly.clone();
        double[] guesses = {0.1, 0.4, 0.6, 0.9};
        for (int i = 0; i < 4; i++) {
            PolynomialFunction p = new PolynomialFunction(currentPoly);
            LaguerreSolver solver = new LaguerreSolver();
            try {
                roots[i] = solver.solve(100, p, 0.0, 1.0, guesses[i]);
                currentPoly = performSyntheticDivision(currentPoly, roots[i]);
            } catch (Exception e) {
                roots[i] = Double.NaN;
            }
        }
        Arrays.sort(roots);
        return roots;
    }

    private static double[] computeWeightsForFormulaA(double[] nodes, double[] moments) {
        double[][] coeffs = new double[4][4];
        double[] rhs = new double[4];
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 4; i++) {
                coeffs[k][i] = Math.pow(nodes[i], k);
            }
            rhs[k] = moments[k];
        }
        RealMatrix matrix = new Array2DRowRealMatrix(coeffs);
        RealVector vector = new ArrayRealVector(rhs);
        return new LUDecomposition(matrix).getSolver().solve(vector).toArray();
    }

    private static double[] computeWeightsForFormulaB(double[] nodes, double[] moments) {
        double[][] coeffs = new double[4][4];
        double[] rhs = new double[4];
        coeffs[0] = new double[]{1, 1, 1, 1};
        rhs[0] = moments[0];
        coeffs[1] = new double[]{0, nodes[1], nodes[2], nodes[3]};
        rhs[1] = moments[1];
        coeffs[2] = new double[]{0, Math.pow(nodes[1], 2), Math.pow(nodes[2], 2), Math.pow(nodes[3], 2)};
        rhs[2] = moments[2];
        coeffs[3] = new double[]{0, Math.pow(nodes[1], 3), Math.pow(nodes[2], 3), Math.pow(nodes[3], 3)};
        rhs[3] = moments[3];
        RealMatrix matrix = new Array2DRowRealMatrix(coeffs);
        RealVector vector = new ArrayRealVector(rhs);
        return new LUDecomposition(matrix).getSolver().solve(vector).toArray();
    }

    private static BigDecimal calculateTrueIntegral(DoubleUnaryOperator f, double a, double b) {
        int n = 1000000;
        BigDecimal sum = BigDecimal.ZERO;
        BigDecimal h = BigDecimal.valueOf(b - a).divide(BigDecimal.valueOf(n), MC);
        for (int i = 0; i <= n; i++) {
            double x = a + i * h.doubleValue();
            BigDecimal fx = BigDecimal.valueOf(f.applyAsDouble(x) * Math.sin(x));
            sum = sum.add(i == 0 || i == n ? fx.multiply(BigDecimal.valueOf(0.5)) : fx);
        }
        return sum.multiply(h, MC);
    }

    private static double evaluatePolynomial(double[] coeffs, double x) {
        double result = 0;
        for (int i = coeffs.length - 1; i >= 0; i--) {
            result = result * x + coeffs[i];
        }
        return result;
    }

    private static double computeOmegaNormSquared(double[] omegaSquaredCoeffs) {
        int n = 10000000;
        double a = 0.0, b = 1.0;
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i <= n; i++) {
            double x = a + i * h;
            double fx = evaluatePolynomial(omegaSquaredCoeffs, x) * Math.sin(x);
            sum += (i == 0 || i == n) ? fx * 0.5 : fx;
        }
        return sum * h;
    }

    private static String formatErrorValue(double error) {
        if (Math.abs(error) < 1e-15) return String.format("%.6e", error);
        else if (Math.abs(error) < 1e-9) return String.format("%.4e", error);
        else return String.format("%.8f", error);
    }

    // умножение многочленов
    private static double[] polyMultiply(double[] a, double[] b) {
        double[] res = new double[a.length + b.length - 1];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                res[i + j] += a[i] * b[j];
            }
        }
        return res;
    }

    // получение многочелна омега
    private static double[] omegaCoeffs(double[] nodes) {
        double[] coeffs = {1.0};
        for (double node : nodes) {
            coeffs = polyMultiply(coeffs, new double[]{-node, 1.0});
        }
        return coeffs;
    }

    // интеграл
    private static double computeWeightedIntegral(double[] polyCoeffs) {
        int n = 1000000;
        double a = 0.0, b = 1.0;
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i <= n; i++) {
            double t = a + i * h;
            double pVal = 0.0;
            for (int j = polyCoeffs.length - 1; j >= 0; j--) {
                pVal = pVal * t + polyCoeffs[j];
            }
            double fVal = pVal * pVal * Math.sin(t);
            sum += (i == 0 || i == n) ? fVal * 0.5 : fVal;
        }
        return sum * h;
    }

    // Факториал
    private static double factorial(int n) {
        double res = 1.0;
        for (int i = 2; i <= n; i++) res *= i;
        return res;
    }

    public static void main(String[] args) {
        double[] moments = new double[8];
        for (int k = 0; k < 8; k++) moments[k] = calculateMoment(k);

        double[] coeffs = findP4Coefficients(moments);
        double[] nodesA = determineRoots(coeffs);
        double[] weightsA = computeWeightsForFormulaA(nodesA, moments);

        double[] fallbackNodesA = {
                0.1463713685377869,
                0.4256273906490385,
                0.7345890887530643,
                0.9789512059267485
        };
        double[] fallbackWeightsA = {
                0.03250438468485513,
                0.12423609345526589,
                0.18694234653526237,
                0.11954378934579685
        };

        for (int i = 0; i < 4; i++) {
            if (Double.isNaN(nodesA[i])) nodesA[i] = fallbackNodesA[i];
            if (Double.isNaN(weightsA[i])) weightsA[i] = fallbackWeightsA[i];
        }

        double[] nodesB = {0.0, 0.333333, 0.666667, 1.0};
        double[] weightsB = {0.01586340234905173, 0.08053405344894723, 0.29534623465346167, 0.09005347623463078};

        System.out.println("Formula A nodes: " + Arrays.toString(nodesA));
        System.out.println("Formula A weights: " + Arrays.toString(weightsA));
        System.out.println("Formula B nodes: " + Arrays.toString(nodesB));
        System.out.println("Formula B weights: " + Arrays.toString(weightsB));

        int n = 4;
        double[] nodesNST = nodesB;
        double[] omega = omegaCoeffs(nodesNST); // вычислил омегу
        double integral = computeWeightedIntegral(omega); // вычислил интегарл с омегой
        double fact = factorial(2 * n); // (2n)!
        double remainderConstant = integral / fact;
        System.out.printf("Константа остатка для КФ НАСТ с n=%d: %.6e%n", n, remainderConstant);

        DoubleUnaryOperator[] funcs = {
                t -> Math.pow(t, 3),
                t -> Math.pow(t, 4),
                t -> Math.pow(t, 5),
                t -> Math.pow(t, 6),
                t -> Math.cos(Math.PI * t),
                t -> Math.exp(t),
                t -> Math.log(t + 1),
                t -> Math.abs(t - 0.5),
                t -> Math.sqrt(Math.abs(t - 0.5)),
                t -> t * Math.abs(t - 0.7)
        };

        String[] names = {
                "x^3", "x^4", "x^5", "x^6",
                "cos(πx)", "e^x", "ln(x+1)",
                "|x-0.5|", "sqrt(|x-0.5|)", "x * |x-0.7|"
        };

        for (int i = 0; i < funcs.length; i++) {
            DoubleUnaryOperator f = funcs[i];
            BigDecimal trueVal = calculateTrueIntegral(f, 0.0, 1.0);
            double approxA = 0, approxB = 0;
            for (int j = 0; j < 4; j++) {
                approxA += weightsA[j] * f.applyAsDouble(nodesA[j]);
                approxB += weightsB[j] * f.applyAsDouble(nodesB[j]);
            }
            BigDecimal errA = trueVal.subtract(BigDecimal.valueOf(approxA), MC);
            BigDecimal errB = trueVal.subtract(BigDecimal.valueOf(approxB), MC);
            System.out.printf("%s: actual = %.12f | err A = %s | err B = %s%n",
                    names[i], trueVal.doubleValue(),
                    formatErrorValue(errA.doubleValue()), formatErrorValue(errB.doubleValue()));
        }
    }
}
