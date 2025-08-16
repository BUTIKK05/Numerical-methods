package org.example.lab3_z5_;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;

public class AdaptiveRK4ConvergenceChart extends Application {

    static double[] f(double x, double[] y) {
        return new double[]{
                y[1],
                y[0] * Math.sin(Math.PI * x) - 2 * y[1]
        };
    }

    static double[] vecAdd(double[] a, double[] b) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) r[i] = a[i] + b[i];
        return r;
    }

    static double[] vecSub(double[] a, double[] b) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) r[i] = a[i] - b[i];
        return r;
    }

    static double[] vecScale(double[] a, double s) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) r[i] = a[i] * s;
        return r;
    }

    static double vecNorm(double[] a) {
        double sum = 0;
        for (double v : a) sum += v * v;
        return Math.sqrt(sum);
    }

    // Шаг RK4
    static double[] rk4Step(double x, double[] y, double h) {
        double[] k1 = f(x, y);
        double[] k2 = f(x + h / 2, vecAdd(y, vecScale(k1, h / 2)));
        double[] k3 = f(x + h / 2, vecAdd(y, vecScale(k2, h / 2)));
        double[] k4 = f(x + h, vecAdd(y, vecScale(k3, h)));

        double[] increment = vecScale(vecAdd(vecAdd(k1, vecScale(k2, 2)), vecAdd(vecScale(k3, 2), k4)), h / 6);
        return vecAdd(y, increment);
    }

    // Шаг RK3
    static double[] rk3Step(double x, double[] y, double h) {
        double[] k1 = f(x, y);
        double[] k2 = f(x + h / 2, vecAdd(y, vecScale(k1, h / 2)));
        double[] k3 = f(x + h, vecAdd(y, vecScale(k2, -h)));

        double[] incr = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            incr[i] = h / 6.0 * (k1[i] + 4 * k2[i] + k3[i]);
        }
        return vecAdd(y, incr);
    }

    static class AdaptiveStepResult {
        double x;
        double[] y;
        double h;
        double localError;

        AdaptiveStepResult(double x, double[] y, double h, double localError) {
            this.x = x;
            this.y = y;
            this.h = h;
            this.localError = localError;
        }
    }

    static AdaptiveStepResult adaptiveRKStep(double x, double[] y, double h, double eps) {
        while (true) {
            double[] yRK4 = rk4Step(x, y, h);
            double[] yRK3 = rk3Step(x, y, h);
            double localError = vecNorm(vecSub(yRK4, yRK3));

            if (localError <= eps) {
                double newH = h;
                if (localError < eps / 4) newH = h * 2;
                return new AdaptiveStepResult(x + h, yRK4, newH, localError);
            } else {
                h = h / 2;
                if (h < 1e-10) {
                    return new AdaptiveStepResult(x, y, h, localError);
                }
            }
        }
    }

    static double exactSolutionU(double x) {
        return Math.exp(-x) * Math.cos(x);
    }

    @Override
    public void start(Stage stage) {
        double x0 = 0.0, xf = 10.0;
        double[] y0 = new double[]{1.0, 0.0};
        double hStart = 0.1;

        List<XYChart.Data<Number, Number>> convergencePoints = new ArrayList<>();

        for (int k = 1; k <= 10; k++) {
            double eps = Math.pow(10, -k);

            double x = x0;
            double[] y = y0.clone();
            double h = hStart;

            while (x < xf) {
                if (x + h > xf) h = xf - x;

                AdaptiveStepResult step = adaptiveRKStep(x, y, h, eps);
                if (step.x == x) break;

                x = step.x;
                y = step.y;
                h = step.h;
            }

            double exact = exactSolutionU(xf);
            double globalError = Math.abs(y[0] - exact);
            convergencePoints.add(new XYChart.Data<>(Math.log10(eps), Math.log10(globalError)));
        }

        NumberAxis xAxis = new NumberAxis();
        xAxis.setLabel("log₁₀(ε) ");

        NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("log₁₀(глобальная ошибка)");

        LineChart<Number, Number> chart = new LineChart<>(xAxis, yAxis);
        chart.setTitle("Сходимость адаптивного RK4 (логарифмическая шкала)");

        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Adaptive RK4");
        series.getData().addAll(convergencePoints);

        chart.getData().add(series);

        VBox root = new VBox(chart);
        Scene scene = new Scene(root, 900, 700);
        stage.setScene(scene);
        stage.setTitle("Логарифмическая диаграмма сходимости");
        stage.show();
    }

    public static void main(String[] args) {
        launch();
    }
}
