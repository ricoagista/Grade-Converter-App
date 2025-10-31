package converter;

public class GradeConverter {

    public Result convert(double nilai) {
        String grade;
        String status;

        if (nilai >= 85 && nilai <= 100) {
            grade = "A";
            status = "Lulus";
        } else if (nilai >= 70) {
            grade = "B";
            status = "Lulus";
        } else if (nilai >= 55) {
            grade = "C";
            status = "Lulus";
        } else if (nilai >= 40) {
            grade = "D";
            status = "Tidak Lulus";
        } else if (nilai >= 0) {
            grade = "E";
            status = "Tidak Lulus";
        } else {
            grade = "-";
            status = "Nilai tidak valid";
        }

        return new Result(grade, status);
    }
}
