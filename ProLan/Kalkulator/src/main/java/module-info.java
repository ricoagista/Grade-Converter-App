module com.kalkulator {
    requires javafx.controls;
    requires javafx.fxml;


    opens com.kalkulator to javafx.fxml;
    exports com.kalkulator;
}