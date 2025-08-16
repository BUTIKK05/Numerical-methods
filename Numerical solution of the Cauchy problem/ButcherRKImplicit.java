package org.example.lab3_z3_z3;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;

public class ButcherRKImplicit extends Application {

    static class State {
        double x, u, v;
        State(double x, double u, double v) {
            this.x = x;
            this.u = u;
            this.v = v;
        }
    }

    // Правая часть: u' = v, v' = u * sin(pi x) - 2v
    static double[] f(double x, double[] y) {
        return new double[]{
                y[1],
                y[0] * Math.sin(Math.PI * x) - 2 * y[1]
        };
    }

    // Неявный метод РК (таблица Бутчера, 2 этапа)
    static double[] butcherRKStep(double x, double[] y, double h) {
        double[][] A = {
                {1.0 / 3, 1.0 / 3},
                {1.0 / 4, 3.0 / 4}
        };
        double[] c = {2.0 / 3, 1.0};
        double[] b = {1.0 / 4, 3.0 / 4};

        int s = 2; // 2 этапа
        double[][] K = new double[s][2];

        // Начальное приближение
        for (int i = 0; i < s; i++) {
            double xi = x + c[i] * h;
            K[i] = f(xi, y);
        }

        // Простые итерации для уточнения (можно заменить Ньютоном)
        for (int iter = 0; iter < 8; iter++) {
            double[][] K_new = new double[s][2];
            for (int i = 0; i < s; i++) {
                double[] y_i = y.clone();
                for (int j = 0; j < s; j++) {
                    for (int k = 0; k < 2; k++) {
                        y_i[k] += h * A[i][j] * K[j][k];
                    }
                }
                K_new[i] = f(x + c[i] * h, y_i);
            }
            K = K_new;
        }

        // Вычисляем y_{n+1}
        double[] yNext = y.clone();
        for (int k = 0; k < 2; k++) {
            for (int i = 0; i < s; i++) {
                yNext[k] += h * b[i] * K[i][k];
            }
        }
        return yNext;
    }

    // Решение всей задачи фиксированным шагом
    static List<State> solveButcherRK(double x0, double[] y0, double xf, double h) {
        List<State> result = new ArrayList<>();
        double x = x0;
        double[] y = y0.clone();

        result.add(new State(x, y[0], y[1]));

        while (x < xf) {
            double[] yNext = butcherRKStep(x, y, h);
            x += h;
            y = yNext;
            result.add(new State(x, y[0], y[1]));
        }
        return result;
    }

    @Override
    public void start(Stage stage) {
        double x0 = 0, xf = 10, h = 0.1;
        double[] y0 = {1.0, 0.0};

        List<State> solution = solveButcherRK(x0, y0, xf, h);
        State last = solution.get(solution.size() - 1);

        NumberAxis xAxis = new NumberAxis(x0, xf, 1);
        xAxis.setLabel("x");

        NumberAxis yAxis = new NumberAxis();
        yAxis.setLabel("u(x), v(x)");

        LineChart<Number, Number> chart = new LineChart<>(xAxis, yAxis);
        chart.setTitle("Неявный метод Рунге–Кутты (таблица Бутчера)");

        XYChart.Series<Number, Number> uSeries = new XYChart.Series<>();
        uSeries.setName("u(x)");

        XYChart.Series<Number, Number> vSeries = new XYChart.Series<>();
        vSeries.setName("v(x)");

        for (State s : solution) {
            uSeries.getData().add(new XYChart.Data<>(s.x, s.u));
            vSeries.getData().add(new XYChart.Data<>(s.x, s.v));
        }

        chart.getData().addAll(uSeries, vSeries);

        Label resultLabel = new Label(String.format("u(10) = %.6f, v(10) = %.6f", last.u, last.v));

        VBox root = new VBox(chart, resultLabel);
        Scene scene = new Scene(root, 900, 600);

        stage.setScene(scene);
        stage.show();
    }

    public static void main(String[] args) {
        launch();
    }
}
