package org.example.lab3_z4_;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;
import javafx.scene.shape.Rectangle;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;

public class RK4VsButcherPoints extends Application {

    static class State {
        double x, u, v;
        State(double x, double u, double v) {
            this.x = x;
            this.u = u;
            this.v = v;
        }
    }

    static double[] f(double x, double u, double v) {
        return new double[]{
                v,
                u * Math.sin(Math.PI * x) - 2 * v
        };
    }

    static double[] f(double x, double[] y) {
        return new double[]{
                y[1],
                y[0] * Math.sin(Math.PI * x) - 2 * y[1]
        };
    }

    static List<State> solveRK4(double x0, double u0, double v0, double xf, double h) {
        List<State> result = new ArrayList<>();
        double x = x0, u = u0, v = v0;
        result.add(new State(x, u, v));
        while (x < xf) {
            double[] k1 = f(x, u, v);
            double[] k2 = f(x + h / 2, u + h / 2 * k1[0], v + h / 2 * k1[1]);
            double[] k3 = f(x + h / 2, u + h / 2 * k2[0], v + h / 2 * k2[1]);
            double[] k4 = f(x + h, u + h * k3[0], v + h * k3[1]);

            u += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
            v += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
            x += h;
            result.add(new State(x, u, v));
        }
        return result;
    }

    static List<State> solveRK4Adaptive(double x0, double u0, double v0, double xf, double h0, double tol) {
        List<State> result = new ArrayList<>();
        double x = x0, u = u0, v = v0;
        double h = h0;

        result.add(new State(x, u, v));

        while (x < xf) {
            if (x + h > xf) h = xf - x;

            double[] k1 = f(x, u, v);
            double[] k2 = f(x + h / 2, u + h / 2 * k1[0], v + h / 2 * k1[1]);
            double[] k3 = f(x + h / 2, u + h / 2 * k2[0], v + h / 2 * k2[1]);
            double[] k4 = f(x + h, u + h * k3[0], v + h * k3[1]);

            double uFull = u + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
            double vFull = v + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

            double h2 = h / 2;

            double[] k1a = f(x, u, v);
            double[] k2a = f(x + h2 / 2, u + h2 / 2 * k1a[0], v + h2 / 2 * k1a[1]);
            double[] k3a = f(x + h2 / 2, u + h2 / 2 * k2a[0], v + h2 / 2 * k2a[1]);
            double[] k4a = f(x + h2, u + h2 * k3a[0], v + h2 * k3a[1]);

            double uHalf = u + h2 / 6 * (k1a[0] + 2 * k2a[0] + 2 * k3a[0] + k4a[0]);
            double vHalf = v + h2 / 6 * (k1a[1] + 2 * k2a[1] + 2 * k3a[1] + k4a[1]);

            double[] k1b = f(x + h2, uHalf, vHalf);
            double[] k2b = f(x + h2 + h2 / 2, uHalf + h2 / 2 * k1b[0], vHalf + h2 / 2 * k1b[1]);
            double[] k3b = f(x + h2 + h2 / 2, uHalf + h2 / 2 * k2b[0], vHalf + h2 / 2 * k2b[1]);
            double[] k4b = f(x + h, uHalf + h2 * k3b[0], vHalf + h2 * k3b[1]);

            double uHalfStep = uHalf + h2 / 6 * (k1b[0] + 2 * k2b[0] + 2 * k3b[0] + k4b[0]);
            double vHalfStep = vHalf + h2 / 6 * (k1b[1] + 2 * k2b[1] + 2 * k3b[1] + k4b[1]);

            double errU = Math.abs(uFull - uHalfStep);
            double errV = Math.abs(vFull - vHalfStep);
            double err = Math.max(errU, errV);

            if (err < tol) {
                x += h;
                u = uHalfStep;
                v = vHalfStep;
                result.add(new State(x, u, v));
            }

            double safety = 0.9;
            if (err == 0) {
                h *= 2;
            } else {
                h *= safety * Math.pow(tol / err, 0.25);
            }

            h = Math.max(h, 1e-6);
            h = Math.min(h, 0.5);
        }

        return result;
    }

    static double[] butcherRKStep(double x, double[] y, double h) {
        double[][] A = {{1.0 / 3, 1.0 / 3}, {1.0 / 4, 3.0 / 4}};
        double[] c = {2.0 / 3, 1.0};
        double[] b = {1.0 / 4, 3.0 / 4};

        int s = 2;
        double[][] K = new double[s][2];
        for (int i = 0; i < s; i++) {
            K[i] = f(x + c[i] * h, y);
        }

        for (int iter = 0; iter < 8; iter++) {
            double[][] Knew = new double[s][2];
            for (int i = 0; i < s; i++) {
                double[] y_i = y.clone();
                for (int j = 0; j < s; j++) {
                    for (int k = 0; k < 2; k++) {
                        y_i[k] += h * A[i][j] * K[j][k];
                    }
                }
                Knew[i] = f(x + c[i] * h, y_i);
            }
            K = Knew;
        }

        double[] yNext = y.clone();
        for (int k = 0; k < 2; k++) {
            for (int i = 0; i < s; i++) {
                yNext[k] += h * b[i] * K[i][k];
            }
        }
        return yNext;
    }

    static List<State> solveButcherRK(double x0, double[] y0, double xf, double h) {
        List<State> result = new ArrayList<>();
        double x = x0;
        double[] y = y0.clone();
        result.add(new State(x, y[0], y[1]));

        while (x < xf) {
            y = butcherRKStep(x, y, h);
            x += h;
            result.add(new State(x, y[0], y[1]));
        }
        return result;
    }

    @Override
    public void start(Stage stage) {
        double x0 = 0, xf = 10, u0 = 1, v0 = 0;

        double initialH = 0.1;
        double tolerance = 1e-6;

        List<State> rk4 = solveRK4Adaptive(x0, u0, v0, xf, initialH, tolerance);

        int N = 200;
        double h = (xf - x0) / N;
        List<State> butcher = solveButcherRK(x0, new double[]{u0, v0}, xf, h);

        State lastRK = rk4.get(rk4.size() - 1);
        State lastB = butcher.get(butcher.size() - 1);

        NumberAxis xAxis = new NumberAxis(x0, xf, 1);
        xAxis.setLabel("x");

        NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("u(x), v(x)");

        LineChart<Number, Number> chart = new LineChart<>(xAxis, yAxis);
        chart.setTitle("Сравнение: RK4 (адаптивный) и метод Бутчера");

        XYChart.Series<Number, Number> rk4SeriesU = new XYChart.Series<>();
        rk4SeriesU.setName("u(x) RK4 (адаптивный)");
        XYChart.Series<Number, Number> rk4SeriesV = new XYChart.Series<>();
        rk4SeriesV.setName("v(x) RK4 (адаптивный)");

        for (State s : rk4) {
            rk4SeriesU.getData().add(new XYChart.Data<>(s.x, s.u));
            rk4SeriesV.getData().add(new XYChart.Data<>(s.x, s.v));
        }

        XYChart.Series<Number, Number> butcherSeriesU = new XYChart.Series<>();
        butcherSeriesU.setName("u(x) Butcher");

        XYChart.Series<Number, Number> butcherSeriesV = new XYChart.Series<>();
        butcherSeriesV.setName("v(x) Butcher");

        for (State s : butcher) {
            butcherSeriesU.getData().add(new XYChart.Data<>(s.x, s.u));
            butcherSeriesV.getData().add(new XYChart.Data<>(s.x, s.v));
        }

        Label label = new Label(String.format("u(xf) RK4 = %.6f, u(xf) Butcher = %.6f", lastRK.u, lastB.u));
        Label label2 = new Label(String.format("v(xf) RK4 = %.6f, v(xf) Butcher = %.6f", lastRK.v, lastB.v));

        Rectangle rect = new Rectangle(20, 20);
        rect.setFill(javafx.scene.paint.Color.RED);

        VBox vbox = new VBox(10, chart, label, label2);
        vbox.setStyle("-fx-padding: 10");

        chart.getData().addAll(rk4SeriesU, rk4SeriesV, butcherSeriesU, butcherSeriesV);

        Scene scene = new Scene(vbox, 800, 600);
        stage.setScene(scene);
        stage.setTitle("RK4 и Butcher Method");
        stage.show();
    }

    public static void main(String[] args) {
        launch();
    }
}
